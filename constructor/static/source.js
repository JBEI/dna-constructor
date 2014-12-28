var node;
var protocolID;
var dotData; // DOT-language representation of the graph

$j(document).ready(function () {
  $j('#exportDiv').hide();
  $j('#primerReport').hide();
  $j('#addFragmentDiv').hide();
  $j('#addFragmentLink').click(function() {
    $j('#fragmentText').val('');
    $j('#addFragmentDiv').toggle();
  });
  $j('#layersText').hide();
  $j('#layersSubmit').hide();
  $j('#layersSubmit').click(function() {
    var layersRange = $j('#layersText').val();
    var layers;
    if (layersRange) {
      if (layersRange.len == 1) {
        layers = parseInt(layersRange);
      }
      
      else {
        layers = layersRange;
      }
      
      window.open('ajax/export/' + protocolID + '?layers=' + layers) 
    }
    $j('#layersText').hide();
    $j('#layersSubmit').hide();
    $j('#instr').hide();
  });
  $j('#instr').hide();
  $j('#nodeInfo').hide();
  $j('#postDiv').hide();
  $j('#loading').hide();
  $j('#exportLink').click(exportClick);
  $j('#nodeClick').click(nodeClick);
  $j('#usePrimerDefaults').click(setPrimerDefaults);
  $j('#useConstructDefaults').click(setConstructDefaults);
  $j('#linkButton').keyup(function(event){
    if (event.keyCode == 13){
      $j('#linkButton').click();
    }
  });
  $j('#sequenceButton').click(onSequenceSend);
  $j('#protocolsLink').click(showProtocols);
  $j('#primerReport').click(generateReport);
});

function setConstructDefaults() {
  $j('#minOverlapText').val(20);
  $j('#maxOverlapText').val(30);
  $j('#maxOligoLengthText').val(80);
  $j('#maxDeviationText').val(10);
}

function setPrimerDefaults() {
  $j('#minLengthText').val(35);
  $j('#maxLengthText').val(45);
  $j('#threeprimeLengthText').val(8);
  $j('#maxThreeprimeText').val(6);
  $j('#numStepsText').val(2);
}

function exportClick() {
  $j('#instr').toggle();
  $j('#layersText').toggle();
  $j('#layersSubmit').toggle();
}

function nodeClick(id) {
  $j('#nodeInfo').jqm({ajax: 'ajax/query/' + id});
  $j('#nodeInfo').jqmShow();
}

function showProtocols() {
  $j('#showProtocols').jqm({ajax: 'ajax/protocols/'});
  $j('#showProtocols').jqmShow();
}

function generateReport() {
  $j.ajax({
    url:'ajax/report/',
    type:'POST',
    data:({id: protocolID, maxPrimerLength: maxLength, minPrimerLength: minLength}),
    success: function(response) {
      respObj = JSON.parse(response);
      $j('#report').html(respObj.stdout);
    }
  });
}

function edgeClick(id) {
  var r = /edgeClick\('(\w+)\s+'\)/g;
  var sendData = dotData.replace(r, "edgeClick('$1')");
  $j.ajax({
    url:'ajax/shrink/',
    type:'POST',
    data:({dotData: sendData, id: id}),
    success: function(response) {
      respObj = JSON.parse(response);
      dotData = respObj.graphData;
      canviz.parse(dotData);
    }
  });
}

function writeResponse(response) {
  respObj = JSON.parse(response);
  
  $j('#loading').hide();
 
  if(respObj.error != false){
    alert(respObj.error);
  } else {
    protocolID = respObj.id;
    dotData = respObj.chartRequest.replace(/\\/g, '');

    canviz = new Canviz('responseDiv');
    canviz.parse(dotData);
    $j('#responseDiv').height($j('[id ^= "canviz"]').height());
    $j('[id ^= "canviz"]').css({'z-index': 1});
    $j('#responseDiv').children('div').css('z-index', 2);
    $j('#responseDiv').show();
    $j('#exportDiv').show();
    $j('#primerReport').show();
  }
}

function writeError(response) {
  $j('#loading').hide();
  $j('#responseDiv').html('Error.');
}

function onSequenceSend() {
  $j('#exportDiv').hide();
  $j('#layersText').hide();
  $j('#layersSubmit').hide();
  $j('#instr').hide();
  $j('#responseDiv').hide();
  sequence = $j('#sequenceText').val().replace(/ /g, '');
  minLength = $j('#minLengthText').val();
  maxLength = $j('#maxLengthText').val();
  threeprimeLength = $j('#threeprimeLengthText').val();
  maxThreeprime = $j('#maxThreeprimeText').val();
  numSteps = $j('#numStepsText').val();
  minOverlap = $j('#minOverlapText').val();
  maxOverlap = $j('#maxOverlapText').val();
  maxOligoLength = $j('#maxOligoLengthText').val();
  maxDeviation = $j('#maxDeviationText').val();
  
  if ($j('#fragmentText').val())
    fragment = $j('#fragmentText').val();
    
  else
    fragment = false;
  
  if (!(minLength && (maxLength && (threeprimeLength && (maxThreeprime && numSteps))))) {
    setPrimerDefaults();
    minLength = $j('#minLengthText').val();
    maxLength = $j('#maxLengthText').val();
    threeprimeLength = $j('#threeprimeLengthText').val();
    maxThreeprime = $j('#maxThreeprimeText').val();
    numSteps = $j('#numStepsText').val();
  }
  
  if (!(minOverlap && (maxOverlap && (maxOligoLength && maxDeviation)))) {
    setConstructDefaults();
    minOverlap = $j('#minOverlapText').val();
    maxOverlap = $j('#maxOverlapText').val();
    maxOligoLength = $j('#maxOligoLengthText').val();
    maxDeviation = $j('#maxDeviationText').val();
  }
  
  $j('#responseDiv').html("");
  
  $j('#loading').show();
  
  $j.ajax({
    url:'ajax/build/',
    type:'POST',
    data:({code: sequence, minLength: minLength, maxLength: maxLength, threeprimeLength: threeprimeLength, 
      maxThreeprime: maxThreeprime, numSteps: numSteps, minOverlap: minOverlap, maxOverlap: maxOverlap, 
      maxOligoLength: maxOligoLength, maxDeviation: maxDeviation, fragment: fragment}),
    success: writeResponse,
    error: writeError
  });
}
