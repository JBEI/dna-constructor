from django.http import HttpResponse, HttpResponseRedirect
from django import forms
from django.utils import simplejson
from django.template import Context, loader
from django.shortcuts import render_to_response
from django.core.servers.basehttp import FileWrapper
from constructor.models import Node, Protocol
from errors import InterpreterError, ConstructError
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from cStringIO import StringIO
import zipfile
import construction
import primers
import interpreter
import os
import re
import pydot
import sys
import subprocess

class SimpleFileForm(forms.Form):
	file = forms.Field(widget=forms.FileInput, required=False)

def index(request):
	return HttpResponseRedirect('/construct/')

def construct(request):
	t = loader.get_template('constructor/construct.html')
	form = SimpleFileForm()
	c = Context({'form': form})
	return HttpResponse(t.render(c))

def gvreturn(request, gv_id):
	f = open(gv_id+'.gv')
	data = f.read()
	response = HttpResponse(data, mimetype='text/plain')
	response['Content-Disposition'] = 'attachment; filename='+gv_id+'.gv'
	return response

def build(request):
	paramsDict = {}
	paramsDict['minPrimerLength'] = (int)(request.POST['minLength']) # get relevant variables from the client's AJAX request
	paramsDict['maxPrimerLength'] = (int)(request.POST['maxLength'])
	paramsDict['threeprimeLength'] = (int)(request.POST['threeprimeLength'])
	paramsDict['maxThreeprime'] = (int)(request.POST['maxThreeprime'])
	paramsDict['numSteps'] = (int)(request.POST['numSteps'])
	paramsDict['minOverlap'] = (int)(request.POST['minOverlap'])
	paramsDict['maxOverlap'] = (int)(request.POST['maxOverlap'])
	paramsDict['maxOligoLength'] = (int)(request.POST['maxOligoLength']) 
	paramsDict['maxDeviation'] = (int)(request.POST['maxDeviation'])
	paramsDict['fragmentList'] = ([request.POST['fragment'].replace(' ','')] if request.POST['fragment'] else False)
	
	try:
		targets, name, id, isExisting = interpreter.parseString(request.POST['code'], paramsDict)

		if isExisting:
			return HttpResponse(simplejson.dumps({'chartRequest':targets.graphString, 'id':targets.id, 'error':False}))
   
		protocol = construction.createProtocol(paramsDict, 'string', targets, name, id)
		reply = {'chartRequest':protocol[0], 'id':protocol[1], 'error':False}
	
		return HttpResponse(simplejson.dumps(reply))

	except InterpreterError as e:
		protocol = Protocol.objects.filter(id=e.protocolId)[0]
		protocol.delete()
		return HttpResponse(simplejson.dumps({'error':e.value}))

	except ConstructError as e:
		protocol = Protocol.objects.filter(id=e.protocolId)[0]
		protocol.delete()
		return HttpResponse(simplejson.dumps({'error':e.value}))

def protocols(request):
	protocols = Protocol.objects.all()
	return render_to_response('constructor/protocols.html', {'protocols':protocols})

def query(request, pk):
	queryNode = Node.objects.get(id=pk)
	return render_to_response('constructor/node.html', {'queryNode':queryNode})

def shrink(request):
	dotData = request.POST['dotData']
	name = request.POST['id']
	graphData = construction.hide_nodes(dotData, name)
	
	reply = {'graphData': graphData}
	return HttpResponse(simplejson.dumps(reply))

def upload(request):
	if request.method == 'POST' and request.FILES['file']:
		destination = open('data/'+request.FILES['file'].name, 'wb+')
		destination.write(request.FILES['file'].read())
		destination.close()
		sequences = construction.parsej5Zip(request.FILES['file'].name)
		return HttpResponse(simplejson.dumps({'sequences': sequences}))
	
	else:
		return HttpResponse(simplejson.dumps({'success': False}))

def report(request):
	protocol = Protocol.objects.get(id=request.POST['id'])
	
	primerList = protocol.getPrimerSequencePairs()
	
	p3process = subprocess.Popen(['/home/primer3-2.2.3/src/primer3_core'], shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	
	p3process.stdin.write('PRIMER_TASK=check_primers\n')
	p3process.stdin.write('PRIMER_MAX_SIZE=' + request.POST['maxPrimerLength'] + '\n')

	for pairDict in primerList:
		p3process.stdin.write('SEQUENCE_TEMPLATE=' + pairDict['target'] + '\n')
		
		if len(pairDict['forwardPrimer']) < 35:
			forwardPrimer = pairDict['forwardPrimer']
		else:
			forwardPrimer = pairDict['forwardPrimer'][len(pairDict['forwardPrimer']) - 35:]
	
		p3process.stdin.write('SEQUENCE_PRIMER=' + forwardPrimer + '\n')
		p3process.stdin.write('=\n')
	
	result, stderr = p3process.communicate()
	return HttpResponse(simplejson.dumps({'stdout':result}))

def export(request, pk):
	requestProtocol = Protocol.objects.get(id=pk)
	
	layers = request.GET['layers'].split(' ')

	if len(layers) > 1:
		format = layers[1]
		layerType = layers[0]
		
	else:
		format = 'fasta'
		layerType = ''.join(layers)
	
	if layerType == 'target':
		nodes = [node for node in Node.objects.filter(protocol=requestProtocol).filter(isRoot=True)]
		for node in nodes:
			if [item.name for item in nodes].count(node.name) > 1:
				for checkNode in nodes[nodes.index(node) + 1:]:
					if checkNode.name == node.name:
						nodes.remove(checkNode)
		
	elif layerType == 'oligos':
		nodes = [node for node in Node.objects.filter(protocol=requestProtocol).filter(isLeaf=True)]
		
		for node in nodes:
			for checkNode in nodes[nodes.index(node) + 1:]:
				if checkNode.name == node.name or checkNode.target == node.target:
					nodes.remove(checkNode)
	
	elif layerType == 1:
		layers = request.GET['layers']
		nodes = [node for node in Node.objects.filter(protocol=requestProtocol).filter(treeDepth=layers)]
		for node in nodes:
			if [item.name for item in nodes].count(node.name) > 1:
				for checkNode in nodes[nodes.index(node) + 1:]:
					if checkNode.name == node.name:
						nodes.remove(checkNode)
	
	else:
		layers = request.GET['layers'].split(' ')[0].split('-')
		layers = [(int)(num) for num in layers]
		layers = range(min(layers), max(layers)+1)
		nodes = [node for node in Node.objects.filter(protocol=requestProtocol).filter(treeDepth__in=layers)]
		
		for node in nodes:
			if [item.name for item in nodes].count(node.name) > 1:
				for checkNode in nodes[nodes.index(node):]:
					if checkNode.name == node.name and [item.name for item in nodes].count(node.name) > 1:
						nodes.remove(checkNode)
						
	if len(request.GET['layers'].split(' ')) > 2: # filter out nodes that aren't a part of the desired target (V1, V2, etc.)
		protocol = request.GET['layers'].split(' ')[2]
		nodesAltered = nodes
		for node in nodes:
			print node
			if node.name.count('Common') == 0 and node.name.count(str(protocol)) == 0:
				nodesAltered.remove(node)
				print len(nodesAltered)
		nodes = nodesAltered
		print nodes, len(nodes)
		
	seqList = []
	primerList = []
	for node in nodes:
		seq = Seq(node.target, generic_dna)
		seqrecord = SeqRecord(seq, id=node.name)
		seqList.append(seqrecord)
		
		if layerType == 'oligos':
			if len(node.leftParentNode.all()) > 0:
				primer = Seq(node.leftParentNode.all()[0].forwardPrimer.sequence, generic_dna)
				record = SeqRecord(primer, id=node.name + '_f_prm')
		
			else:
				primer = Seq(node.rightParentNode.all()[0].reversePrimer.sequence, generic_dna)
				record = SeqRecord(primer, id=node.name + '_r_prm')
				
			primerList.append(record)
		
		elif not(node.isLeaf):
			fprimer = Seq(node.forwardPrimer.sequence, generic_dna)
			frecord = SeqRecord(fprimer, id=node.name+'_Primer_Forward')
			rprimer = Seq(node.reversePrimer.sequence, generic_dna)
			rrecord = SeqRecord(rprimer, id=node.name+'_Primer_Reverse')
			primerList.append(frecord)
			primerList.append(rrecord)
	
	if format == 'plain':
		seqFilename = 'data/' + str(requestProtocol.name) + '_sequences.txt'
		prmFilename = 'data/' + str(requestProtocol.name) + '_primers.txt'
		seqPath = os.path.join(os.getcwd(), seqFilename)
		prmPath = os.path.join(os.getcwd(), prmFilename)
		seqOutputHandle = open(seqPath, 'w')
		prmOutputHandle = open(prmPath, 'w')
		
		for rec in seqList:
			seqOutputHandle.write(rec.id + '\t' + str(rec.seq) + '\t' + str(len(str(rec.seq))) + '\n')
		
		seqOutputHandle.close()
		
		for rec in primerList:
			prmOutputHandle.write(rec.id + '\t' + str(rec.seq) + '\t' + str(len(str(rec.seq))) + '\n')
		
		prmOutputHandle.close()
	
	else:
		seqFilename = 'data/' + str(requestProtocol.name) + '_sequences.fasta'
		prmFilename = 'data/' + str(requestProtocol.name) + '_primers.fasta'
		seqPath = os.path.join(os.getcwd(), seqFilename)
		prmPath = os.path.join(os.getcwd(), prmFilename)
		seqOutputHandle = open(seqPath, 'w')
		prmOutputHandle = open(prmPath, 'w')
		
		SeqIO.write(seqList, seqOutputHandle, 'fasta')
		SeqIO.write(primerList, prmOutputHandle, 'fasta')
		seqOutputHandle.close()
		prmOutputHandle.close()
	
	zipName = 'data/' + str(requestProtocol.name) + '_' + layerType + '.zip'
	zipName = os.path.join(os.getcwd(), zipName)
	zip = zipfile.ZipFile(zipName, mode='w')
	zip.write(seqPath, seqFilename)
	zip.write(prmPath, prmFilename)
	zip.close()
								  
	wrapper = FileWrapper(file(zipName))
	response = HttpResponse(wrapper, content_type='application/octet-stream')
	response['Content-Length'] = os.path.getsize(zipName)
	response['Content-Disposition'] = 'attachment; filename="protocol_' + zipName[6:] + '"'
	return response
