// Auto-size the embedding iframe to this artifact's natural height.
// Parent listens for { type: 'artifact-resize', height } and sets iframe.style.height.
(function () {
  if (window.parent === window) return;

  let lastHeight = 0;

  function measure() {
    const doc = document.documentElement;
    const body = document.body;
    return Math.max(
      body.scrollHeight, body.offsetHeight,
      doc.clientHeight, doc.scrollHeight, doc.offsetHeight
    );
  }

  function post() {
    const h = measure();
    if (h === lastHeight) return;
    lastHeight = h;
    window.parent.postMessage({ type: 'artifact-resize', height: h }, '*');
  }

  if (typeof ResizeObserver !== 'undefined') {
    const ro = new ResizeObserver(post);
    ro.observe(document.documentElement);
    ro.observe(document.body);
  } else {
    setInterval(post, 250);
  }

  window.addEventListener('load', post);
  window.addEventListener('resize', post);
  document.addEventListener('DOMContentLoaded', post);
})();
