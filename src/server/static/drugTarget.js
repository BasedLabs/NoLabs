var spinnerLoadingExperimentsEnable = function () {
    $('#submitInference').attr('disabled', true);
    $('#spinnerExperiments').removeClass('invisible');
}

var spinnerLoadingExperimentsDisable = function () {
    $('#submitInference').attr('disabled', false);
    $('#spinnerExperiments').addClass('invisible');
}

var loadExistingExperiments = function () {
    $.ajax({
        url: 'load-experiments',
        type: 'GET',
        async: true,
        cache: false,
        contentType: false,
        processData: false,
        success: function (experiments) {

        }
    });
}