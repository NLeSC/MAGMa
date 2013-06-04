Ext.define('Esc.magmaweb.view.Viewport', {
  extend: 'Ext.container.Viewport',
  layout: 'border',
  requires: [
    'Ext.panel.Panel',
    'Ext.layout.container.Border',
    'Ext.Img',
    'Ext.toolbar.Spacer',
    'Ext.container.ButtonGroup',
    'Esc.magmaweb.view.fragment.Tree',
    'Esc.magmaweb.view.metabolite.Panel',
    'Esc.magmaweb.view.scan.Panel'
  ],
  initComponent: function() {
      this.items = [{
          region: 'north',
          layout: {
            type: 'hbox',
            align: 'middle',
            padding: 2
          },
          items: [{
              xtype: 'buttongroup',
              columns: 3,
              items: [{
                  text: 'Home',
                  href: '../',
                  hrefTarget: '_self'
              }, {
                  text: 'Help',
                  href: '../help'
              }, {
                  text: 'Annotate',
                  tooltip: 'Annotate all structures',
                  id: 'annotateaction',
                  iconCls: 'icon-annot',
                  disabled: true
              }, {
                  text: 'Workspace',
                  tooltip: 'My settings and jobs',
                  href: '../workspace',
                  hrefTarget: '_self'
              }, {
                  text: 'Logout',
                  href: '../logout',
                  hrefTarget: '_self'
              }, {
                  text: 'Login',
                  href: "../login",
                  hrefTarget: '_self'
              }, {
                  text: 'Information',
                  tooltip: 'Information about input parameters',
                  action: 'information'
              }]
          }, {
            xtype:'tbspacer',
            flex:1 // aligns buttongroup right
          }, {
            xtype: 'component',
            cls: 'x-title',
            html: '<a href=".." data-qtip="<b>M</b>s <b>A</b>nnotation based on in silico <b>G</b>enerated <b>M</b>et<b>a</b>bolites">MAGMa</a>'
          }, {
            xtype:'tbspacer',
            flex:1 // aligns buttongroup right
          }, {
            xtype: 'component',
            contentEl: 'logos'
          }]
      }, {
          // master side
          region: 'center',
          layout: 'border',
          border: false,
          items:[{
            region:'center',
            border: false,
            xtype: 'metabolitepanel'
          },{
            region:'south',
            hideCollapseTool: true,
            collapsible: true,
            height: '50%',
            split: true,
            xtype: 'scanpanel',
            border: false
          }]
      }, {
          // detail side
          region: 'east',
          split: true,
          collapsible: true,
          layout: 'border',
          width: 600,
          hideCollapseTool: true,
          border: false,
          preventHeader: true,
          items:[{
            region: 'center',
            xtype: 'fragmenttree',
            border: false
          }, {
            region:'south',
            height: '50%',
            split: true,
            collapsible: true,
            hideCollapseTool: true,
            border: false,
            title: 'Scans',
            id: 'mspectrapanel',
            layout: {
                type: 'vbox',
                align: 'stretch'
            },
            defaults: {
                flex: 1,
                layout:'fit',
                border: false
            }
          }]
      }];
      this.callParent();
  }
});
