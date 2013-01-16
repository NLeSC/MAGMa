Installation
============

Each sub-project has it own installation instructions.

Import into Eclipse
-------------------

To get the subprojects as seperate projects into Eclipse do:

1. Import 'Projects from Git'
2. Use URI with 'git@github.com:NLeSC/MAGMa.git'
3. Import all branches
4. Use default destination '~/git/MAGMa'
5. Import 'Projects from Git'
6. Use Local with '~/git/MAGMa'
7. Select import as general project
8. Open working directory tree and select 'emetabolomics_site'
9. Repeat steps 5..8 for 'job' and 'web', rename projects by prefixing 'magma'
10. On 'job' and 'web' project folder right click and select 'PyDev>Set as PyDev Project'
11. On 'job' and 'web' project folder right click and select 'Properties>PyDev PYTHONPATH'
12. Press 'Add source folder' and select project folder
13. Install virtualenv python with packages required by job and web and run 'python setup.py develop' and setup python interpreter.
14. Import Maven>'Existing Maven Projects'
15. Set Root Directory to ~/git/MAGMa/joblauncher
16. On 'joblauncher' project folder right click and select 'Team>Share project'
17. Select git select magma git repo
