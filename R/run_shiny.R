#' @export
mgmwm_gui = function(){
  appDir = system.file("mgmwmapp", package = "demo")
  shiny::runApp(appDir, display.mode = "normal")
}
