// functions for making a user's choices persist across page refreshes
function makePersistent(item, id) {
  // get the stored value and the element
  var storedValue = sessionStorage.getItem(item);
  var element = document.getElementById(id);

  // if a stored value exists, set the element's value
  if (storedValue !== null) {
    element.value = storedValue;
  }

  // add event listener to persist changes
  element.addEventListener("change", function () {
    var value = element.value;
    sessionStorage.setItem(item, value);
  });
}

function makePersistentCheckbox(item, id) {
  // set alignSmarts checkbox if it was previously set
  var storedValue = sessionStorage.getItem(item);
  var element = document.getElementById(id);

  if (storedValue !== null) {
    element.checked = storedValue === "true";
  }

  // add event listener to persist checkbox
  element.addEventListener("change", function () {
    sessionStorage.setItem(item, element.checked);
  });
}

function makePersistentSelect(item, id) {
  var defOption = id + "Default";
  // Get the select element and the default selected option
  var element = document.getElementById(id);
  var defaultOption = document.getElementById(defOption);

  // Check if there is a saved value in storage
  var savedSize = sessionStorage.getItem(item);

  // Set the selected option based on the saved value or default
  defaultOption.selected = true; // Ensure one is always selected
  if (savedSize) {
    // Find the option with the saved value and select it
    var selectedOption = Array.from(element.options).find(
      (option) => option.value === savedSize
    );
    if (selectedOption) {
      selectedOption.selected = true;
    }
  }
  // Add an event listener to save the selected option when it changes
  element.addEventListener("change", function () {
    sessionStorage.setItem(item, element.value);
  });
}

// setup listeners
function setListeners() {
  var modal = document.getElementById("myModal");
  var modalImg = document.getElementById("img01");
  const captionText = document.getElementById("caption");
  var images = document.querySelectorAll(".molimg");
  for (var i = 0; i < images.length; i++) {
    images[i].onclick = function () {
      modal.style.display = "block";
      modalImg.src = "/" + this.getAttribute("value");
      modalImg.alt = this.getAttribute("name");
      captionText.innerHTML = this.getAttribute("name");
    };
  }

  modal.onclick = function () {
    modalImg.className += " out";
    setTimeout(function () {
      modal.style.display = "none";
      modalImg.className = "modal-content";
    }, 400);
  };

  window.onkeyup = function (event) {
    if (event.keyCode == 27) {
      modalImg.className += " out";
      setTimeout(function () {
        modal.style.display = "none";
        modalImg.className = "modal-content";
      }, 400);
    }
  };

  document.getElementById("reset-button").onclick = function () {
    reset();
  };

  document.getElementById("help-button").onclick = function () {
    help();
  };
}

// show/hide input options based on selected format
/* there is some flickering when the page loads, but to fix 
  we'd need to use a server-side solution to render the page before sending it to the client*/
function showHideOptions() {
  var molfmt = document.getElementById("molFmt");
  var options = document.getElementById("user-options-csv");

  function checkMolfmtValue() {
    if (molfmt.value === "SMILES") {
      options.style.display = "block";
      sessionStorage.setItem("optionsDisplay", "block");
    } else {
      options.style.display = "none";
      sessionStorage.setItem("optionsDisplay", "none");
    }
  }

  // Run checkMolfmtValue when the page loads
  if (sessionStorage.getItem("optionsDisplay")) {
    options.style.display = sessionStorage.getItem("optionsDisplay");
  } else {
    checkMolfmtValue();
  }

  // Run checkMolfmtValue when the value of molfmt changes
  molfmt.addEventListener("change", checkMolfmtValue);
}

// // Optional, update the width when the dropdown selection changes
// document.getElementById("mols-row").addEventListener("change", updateImageWidths);

