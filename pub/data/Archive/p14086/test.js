// Create variables for code input field & submit button
var codeInput = document.getElementById('codeField');
var submitBtn = document.getElementById('codeSubmit');
// Hide error message
$('.error').hide();

function submitCode(code) {
  var code = codeInput.value;
  // Search for a match between regular expression and user input
  var regEx = /^[0-9]+$/.test(code);
  console.log(regEx);
  // If code length is not 10 or regular expression doesn't match
  if(code.length != 10 || regEx === false) {
    console.log("Code not accepted");
    // Show error message
    $('.error').show();
  } else {
      console.log("Code accepted");
      // Hides error message if code was previously input incorrectly
      $('.error').hide();
      //window.location="https://www.google.com";
    }
  }

$('#codeSubmit').click( function() {
  submitCode();
})
