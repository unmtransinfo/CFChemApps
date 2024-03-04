// adjust shown images
function adjustImageGrid() {
  // this function is used to center the image grid if there's only one row of images
  // Get the image grid element
  var imageGrid = document.getElementById("imageGrid");

  // Get the height of the image grid and the height of a single image
  var gridHeight = imageGrid.offsetHeight;
  var querySelector = imageGrid.querySelector(".image-content");
  if (querySelector === null) {
    // no images to display, don't need to do anything more
    return;
  }
  var imageHeight = querySelector.offsetHeight;

  // If the height of the image grid is less than or equal to the height of a single image (plus margins), it means there's only one row
  if (gridHeight <= imageHeight + 2) {
    // assuming 1px margin
    imageGrid.style.justifyContent = "center";
  } else {
    imageGrid.style.justifyContent = "flex-start";
  }
}
function adjustImageContentWidths() {
  // adjust the width of image-content to match the width of the images
  var imageContents = document.getElementsByClassName("image-content");
  for (var i = 0; i < imageContents.length; i++) {
    var image = imageContents[i].getElementsByTagName("img")[0];
    var padding = 12;
    var width = image.offsetWidth + padding;
    imageContents[i].style.maxWidth = width + "px";
  }
}
