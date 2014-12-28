from django.db import models

# Class for a node in a reaction tree.
class Node(models.Model):
    name = models.CharField(max_length=1000)
    isRoot = models.BooleanField(default=False)
    isLeaf = models.BooleanField(default=False)
    isNaturalFragment = models.BooleanField(default=False)
    isCommon = models.BooleanField(default=False)
    treeDepth = models.IntegerField()
    protocol = models.ForeignKey('Protocol', related_name='Nodes', null=True, blank=True)    
    target = models.CharField(max_length=10000, null=True, blank=True)
    forwardPrimer = models.ForeignKey('Primer', related_name='forward_owners', null=True, blank=True)
    reversePrimer = models.ForeignKey('Primer', related_name='reverse_owners', null=True, blank=True)
    startIndex = models.IntegerField()
    divPoint = models.IntegerField(null=True, blank=True)
    leftReaction = models.ForeignKey('Node', related_name='leftParentNode', null=True, blank=True)
    rightReaction = models.ForeignKey('Node', related_name='rightParentNode', null=True, blank=True)
    
    def protocolID(self):
        return str(self.protocol.id)
    
    def __unicode__(self):
        return ('(' + self.protocolID() + ') ' + self.name)

    def getAllChildPrimerSequencePairs(self, primers):
        if self.leftReaction:
            primers = self.leftReaction.getAllChildPrimerSequencePairs(primers)
        if self.rightReaction:
            primers = self.rightReaction.getAllChildPrimerSequencePairs(primers)

        if self.forwardPrimer:
            appendDict = {'target':self.target, 'forwardPrimer':self.forwardPrimer.sequence, 'reversePrimer':self.reversePrimer.sequence}
            primers.extend([appendDict])
            return primers
        
        else:
            return primers

class Primer(models.Model):
    sequence = models.CharField(max_length=500)
    size = models.IntegerField()
    
    def protocolID(self):
        if len(self.forward_owners.all()) > 0:
            return self.forward_owners.all()[0].protocolID()
        
        elif len(self.reverse_owners.all()) > 0:
            return self.reverse_owners.all()[0].protocolID()
        
        else:
            return '0'
    
    def __unicode__(self):
        return self.sequence
    
class Protocol(models.Model):
    name = models.CharField(max_length=1000, null=True, blank=True)
    target = models.CharField(max_length=10000, null=True, blank=True)
    definition = models.CharField(max_length=10000, null=True, blank=True)
    params = models.CharField(max_length=1000)
    rootNode = models.ForeignKey('Node', related_name='Root Node', null=True, blank=True)
    graphString = models.CharField(max_length=100000, null=True, blank=True)
    timeCreated = models.DateTimeField(auto_now_add=True)
    
    def rootSequence(self):
        string = self.rootNode.target[0:20] + ' ... ' + self.rootNode.target[len(self.rootNode.target) - 20:]
        return string
    
    def __unicode__(self):
        return ('Protocol ' + str(self.id))

    def getPrimerSequencePairs(self):
        return self.rootNode.getAllChildPrimerSequencePairs([])
    
class SequenceVariable(models.Model):
    """Class for user-defined sequences."""
    name = models.CharField(max_length=1000, null=True, blank=True)
    sequence = models.CharField(max_length=10000, null=True, blank=True)
    protocol = models.ForeignKey('Protocol', related_name='sequenceOwner')
