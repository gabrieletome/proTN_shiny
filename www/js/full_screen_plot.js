/*protn*/
function showFullscreenPlot(id) {
  const container = document.getElementById('fullscreen-container');
  container.classList.add('show');
  Shiny.setInputValue('fullscreen_trigger', id, {priority: 'event'});
}
function hideFullscreenPlot() {
  const container = document.getElementById('fullscreen-container');
  container.classList.remove('show');
  Shiny.setInputValue('fullscreen_trigger', null, {priority: 'event'});
}

/*phosprotn*/
function showFullscreenPlot_phos(id) {
  const container = document.getElementById('fullscreen-container-phos');
  container.classList.add('show');
  Shiny.setInputValue('fullscreen_trigger_phos', id, {priority: 'event'});
}
function hideFullscreenPlot_phos() {
  const container = document.getElementById('fullscreen-container-phos');
  container.classList.remove('show');
  Shiny.setInputValue('fullscreen_trigger_phos', null, {priority: 'event'});
}

/*phosprotn with proteome*/
function showFullscreenPlot_phos_protn(id) {
  const container = document.getElementById('fullscreen-container-phos-protn');
  container.classList.add('show');
  Shiny.setInputValue('fullscreen_trigger_phos_protn', id, {priority: 'event'});
}
function hideFullscreenPlot_phos_protn() {
  const container = document.getElementById('fullscreen-container-phos-protn');
  container.classList.remove('show');
  Shiny.setInputValue('fullscreen_trigger_phos_protn', null, {priority: 'event'});
}

/*interactn*/
function showFullscreenPlot_interactn(id) {
  const container = document.getElementById('fullscreen-container-interactn');
  container.classList.add('show');
  Shiny.setInputValue('fullscreen_trigger_interactn', id, {priority: 'event'});
}
function hideFullscreenPlot_interactn() {
  const container = document.getElementById('fullscreen-container-interactn');
  container.classList.remove('show');
  Shiny.setInputValue('fullscreen_trigger_interactn', null, {priority: 'event'});
}