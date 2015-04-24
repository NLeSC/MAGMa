/**
 * Ajax store which aborts all pending requests when a new request is done.
 */
Ext.define('Esc.magmaweb.store.proxy.AjaxSingle', {
  extend: 'Ext.data.proxy.Ajax',
  alias: 'proxy.ajaxsingle',
  doRequest: function(operation, callback, scope) {
    var writer = this.getWriter(),
      request = this.buildRequest(operation);

    if (operation.allowWrite()) {
      request = writer.write(request);
    }

    Ext.apply(request, {
      binary: this.binary,
      headers: this.headers,
      timeout: this.timeout,
      scope: this,
      callback: this.createRequestCallback(request, operation, callback, scope),
      method: this.getMethod(request),
      disableCaching: false // explicitly set it to false, ServerProxy handles caching
    });

    // cancel any running requests
    if (this.activeRequest) {
      Ext.Ajax.abort(this.activeRequest);
    }

    this.activeRequest = Ext.Ajax.request(request);

    return request;
  }
});