// file read functions
function getExtension(filename) {
  var parts = filename.split(".");
  return parts[parts.length - 1];
}

function readFile() {
  document.getElementById("infile").addEventListener("change", function () {
    console.log("hi");
    // Check if a file was selected
    var file = this.files[0];
    if (file) {
      // Check if the file size exceeds the limit (5 MB)
      var fileSizeMB = file.size / (1024 * 1024);
      if (fileSizeMB > 5) {
        alert(
          "The file size exceeds the limit (5 MB). Please select a smaller file."
        );
        reset((resetOptions = false));
        return;
      }
      // check file extension
      var extension = "." + getExtension(file.name);
      var acceptedFileTypes = this.accept.split(",").map((item) => item.trim());
      console.log(acceptedFileTypes);
      if (!acceptedFileTypes.includes(extension)) {
        var msg = `The detected file type (${extension}) is not supported.`;
        alert(msg);
        reset((resetOptions = false));
        return;
      }

      // Read the file
      var file_reader = new FileReader();
      file_reader.onload = function () {
        document.getElementById("textarea").textContent = file_reader.result;
      };
      file_reader.readAsText(file);
    }
  });
}

function checkPastedChars() {
  // don't allow for user to paste more than maxLength chars
  var textarea = document.getElementById("textarea");
  var maxLength = 5242880; // 5 MB worth of characters

  textarea.addEventListener("input", function () {
    if (this.value.length > maxLength) {
      alert(
        "The input exceeds the limit (" +
          maxLength +
          " characters). It will be truncated."
      );
      this.value = this.value.substring(0, maxLength);
    }
  });
}
