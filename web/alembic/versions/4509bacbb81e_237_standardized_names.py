"""237-standardized-names

Revision ID: 4509bacbb81e
Revises: 185259a481ee
Create Date: 2015-03-20 15:19:02.909799

"""
from alembic import op

# revision identifiers, used by Alembic.
revision = '4509bacbb81e'
down_revision = '185259a481ee'


def upgrade():
    """
    Sqlite can only rename table or add columns.
    As described at http://sqlite.org/lang_altertable.html
    To rename metid 2 molid tables must be replaced.

    The new column smiles is empty.
    """
    connection = op.get_bind().connection
    sql = """
        CREATE TABLE molecules (
            molid INTEGER NOT NULL,
            mol VARCHAR,
            reactionsequence VARCHAR,
            inchikey14 VARCHAR,
            smiles VARCHAR,
            formula VARCHAR,
            predicted BOOLEAN,
            name VARCHAR,
            nhits INTEGER,
            mim FLOAT,
            natoms INTEGER,
            logp FLOAT,
            refscore FLOAT,
            reference VARCHAR,
            PRIMARY KEY (molid),
            UNIQUE (inchikey14),
            CHECK (predicted IN (0, 1))
        );

        INSERT INTO molecules (
            molid, mol, reactionsequence, formula, predicted, name,
            nhits, mim, natoms, logp, refscore, reference, inchikey14
        ) SELECT
            metid, mol, reactionsequence,
            molformula,
            NOT isquery,
            origin,
            nhits, mim, natoms, logp,
            probability,
            reference,
            smiles
        FROM metabolites;

        ALTER TABLE peaks RENAME TO old_peaks;

        DROP INDEX ix_peaks_assigned_metid;

        CREATE TABLE peaks (
            scanid INTEGER NOT NULL,
            mz FLOAT NOT NULL,
            intensity FLOAT,
            assigned_molid INTEGER,
            PRIMARY KEY (scanid, mz),
            FOREIGN KEY(scanid) REFERENCES scans (scanid),
            FOREIGN KEY(assigned_molid) REFERENCES molecules (molid)
        );

        INSERT INTO peaks SELECT * FROM old_peaks;

        CREATE INDEX ix_peaks_assigned_molid ON peaks (assigned_molid);

        ALTER TABLE reactions RENAME TO reactions_old;

        CREATE TABLE reactions (
            reactid INTEGER NOT NULL,
            reactant INTEGER,
            product INTEGER,
            name VARCHAR,
            PRIMARY KEY (reactid),
            FOREIGN KEY(reactant) REFERENCES molecules (molid),
            FOREIGN KEY(product) REFERENCES molecules (molid)
        );
        INSERT INTO reactions SELECT * FROM reactions_old;


        ALTER TABLE fragments RENAME TO fragments_old;

        DROP INDEX ix_fragments_metid;
        DROP INDEX ix_fragments_scanid;

        CREATE TABLE fragments (
            fragid INTEGER NOT NULL,
            molid INTEGER,
            scanid INTEGER,
            mz FLOAT,
            mass FLOAT,
            score FLOAT,
            parentfragid INTEGER,
            atoms VARCHAR,
            deltah FLOAT,
            deltappm FLOAT,
            smiles VARCHAR,
            formula VARCHAR,
            PRIMARY KEY (fragid),
            FOREIGN KEY(scanid, mz) REFERENCES peaks (scanid, mz),
            FOREIGN KEY(molid) REFERENCES molecules (molid),
            FOREIGN KEY(scanid) REFERENCES scans (scanid),
            FOREIGN KEY(parentfragid) REFERENCES fragments (fragid)
        );

        INSERT INTO fragments SELECT * FROM fragments_old;

        CREATE INDEX ix_fragments_molid ON fragments (molid);
        CREATE INDEX ix_fragments_scanid ON fragments (scanid);

        DROP TABLE reactions_old;
        DROP TABLE fragments_old;
        DROP TABLE old_peaks;
        DROP TABLE metabolites;
    """
    connection.executescript(sql)


def downgrade():
    """Can never go back"""
    pass
