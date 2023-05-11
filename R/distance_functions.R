deg2rad <- function(deg){
  return(deg * pi / 180)
}

rad2deg <- function(rad){
  return(rad * 180 / pi)
}

distanceEarth <- function(lat1d, lon1d, lat2d, lon2d, radius) {
  lat1r = deg2rad(lat1d)
  lon1r = deg2rad(lon1d)
  lat2r = deg2rad(lat2d)
  lon2r = deg2rad(lon2d)
  u = sin((lat2r - lat1r)/2);
  v = sin((lon2r - lon1r)/2);
  return(2.0 * radius * sqrt(u^2 + cos(lat1r) * cos(lat2r) * v^2))
}

calculateDistance <- function(x1, y1, x2, y2, radius) {
  if(is.null(radius)){
    return(sqrt((x2 - x1)^2 + (y2 - y1)^2))
  }else{
    return (distanceEarth(x1, y1, x2, y2, radius))
  }
}
