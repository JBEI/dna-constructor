from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import pairwise2
import logging
import time
import re
import sys

def main(paramsDict, seq_type='string', sequence='none', forwardOnly=False, reverseOnly=False):
    """This function reads parameters from a file and finds the optimum forward and reverse primer sizes."""
    
    template = parse_sequence(sequence, seq_type)

    paramsDict['seq_one'] = template
    paramsDict['maxThreeprimeGap'] = 0
    
    if not(reverseOnly):
        bestForwardSize, bestForwardScore, bestForwardValid, forwardData = get_best_primer(
                            paramsDict['seq_one'], paramsDict['minPrimerLength'], paramsDict['maxPrimerLength'],
                            paramsDict['threeprimeLength'], paramsDict['maxThreeprime'], paramsDict['maxThreeprimeGap'],
                            paramsDict['numSteps'])
    
    if not(forwardOnly):
        paramsDict['seq_one'] = template.reverse_complement()
        
        bestReverseSize, bestReverseScore, bestReverseValid, reverseData = get_best_primer(
                            paramsDict['seq_one'], paramsDict['minPrimerLength'], paramsDict['maxPrimerLength'],
                            paramsDict['threeprimeLength'], paramsDict['maxThreeprime'], paramsDict['maxThreeprimeGap'],
                            paramsDict['numSteps'])
    
    if not(forwardOnly) and not(reverseOnly):
        if bestForwardValid and bestReverseValid:
            bestPairSizes, bestPairScores = check_for_dimers(template, forwardData, reverseData)
            bestForwardSize = bestPairSizes[0]
            bestReverseSize = bestPairSizes[1]
            bestForwardScore = bestPairScores[0]
            bestReverseScore = bestPairScores[1]
            
            #print('optimal forward primer size is ' + str(bestForwardSize) + ', with a score of ' +str(bestForwardScore)+ '\n')
            #print('optimal reverse primer size is ' + str(bestReverseSize) + ', with a score of ' +str(bestReverseScore))
            
            return (bestForwardSize, bestForwardScore, bestReverseSize, bestReverseScore)
            
        elif bestForwardValid and not(bestReverseValid):
            #print('optimal forward primer size is ' + str(bestForwardSize) + ', with a score of ' +str(bestForwardScore)+ '\n')
            print('no valid reverse primer size in given size range')
            return (bestForwardSize, bestForwardScore, False, False)
        
        elif not(bestForwardValid) and bestReverseValid:
            print('no valid reverse primer size in given size range')
            #print('optimal reverse primer size is ' + str(bestReverseSize) + ', with a score of ' +str(bestReverseScore))
            return (False, False, bestReverseSize, bestReverseScore)
            
        else:
            print('no valid forward or reverse primer sizes in given size range')
            return (False, False, False, False)
        
    if forwardOnly:
        if bestForwardValid:
            return (bestForwardSize, bestForwardScore)
        
        else:
            return (False, False)
        
    if reverseOnly:
        print 'bestReverseValid: ' + str(bestReverseValid)
        if bestReverseValid:
            return (bestReverseSize, bestReverseScore)
        
        else:
            return (False, False)
    
def get_best_primer(seq_one, minPrimerLength, maxPrimerLength, threeprimeLength, maxThreeprime, maxThreeprimeGap, numSteps):
    """Function to find the best forward primer for a given template sequence."""
    primerSize = minPrimerLength
    seq_two = Seq((str(seq_one)[0:primerSize]).upper(), generic_dna)
    
    score = 0
    valid = True
    data = {'length':[], 'score':[], 'valid':[], 'ratio':[]}
    
    subSeqs = get_sub_sequences(numSteps, seq_one[len(seq_two):], seq_two)
    
    for subSeq in subSeqs:
        
        (subScore, subValid) = find_primer_score(numSteps, threeprimeLength, maxThreeprime, maxThreeprimeGap, subSeq, seq_two)
        score += subScore
        
        if subValid == False:
            valid = False
            break
        
    ratio = (float)(score) / (float)(primerSize)
    
    if valid:
        data = {'length':[primerSize], 'score':[score], 'valid':[valid], 'ratio':[ratio]} # store primer lengths, scores, and whether or not each length is valid in a dict
    
    bestPrimerSize = primerSize
    bestPrimerScore = score
    bestPrimerValid = valid
    
    for i in range(minPrimerLength + 1, maxPrimerLength + 1):
        seq_two = Seq((str(seq_one)[0:i]).upper(), generic_dna)
        subSeqs = get_sub_sequences(numSteps, seq_one[len(seq_two):], seq_two)
        
        score = 0
        valid = True
        
        for subSeq in subSeqs:
        
            #print(subSeq + ':\n')
            
            (subScore, subValid) = find_primer_score(numSteps, threeprimeLength, maxThreeprime, maxThreeprimeGap, subSeq, seq_two)
            score += subScore
            
            if subValid == False:
                valid = False
                break
            
        ratio = (float)(score) / (float)(i)
        
        if valid:
            data['length'].append(i)
            data['score'].append(score)
            data['valid'].append(valid)
            data['ratio'].append(ratio)
        
        if (score < bestPrimerScore or not(bestPrimerValid)) and valid: # save current size as best option if we have no valid primers or it is the lowest score yet
            bestPrimerSize = i
            bestPrimerScore = score
            bestPrimerValid = True
        
        #print('with primer size ' +str(i)+ ', valid = ' +str(valid)+ ', score = ' +str(score)+ ', ratio = ' +str(ratio)+'\n')
    
    return bestPrimerSize, bestPrimerScore, bestPrimerValid, data
    
def find_primer_score(numSteps, threeprimeLength, maxThreeprime, maxThreeprimeGap, seq_one, seq_two):
    """Returns alignment score and validity of a given primer/template combination."""
    totalScore = 0.0
    
    locationsStr = []
    
    locationScore = 0.0
    mostThreeprimeMatchesSoFar = 0
    alignment = pairwise2.align.localms(seq_two, seq_one, 2, -1, -.5, -.1, one_alignment_only=True)
    
    valid = True
    
    if len(seq_one) <= len(seq_two):
        valid = False
        return(0, valid)
        
    for align1, align2, score, begin, end in alignment:
        
        matchLocations, primerLocations = get_indices_from_alignment(align1, align2)
        
        threeprimeStart = primerLocations[len(primerLocations) - threeprimeLength]
        
        threeprimeMatches = get_threeprime_matches(matchLocations, threeprimeStart, threeprimeLength, maxThreeprimeGap)
                    
        if threeprimeMatches > maxThreeprime: # if we have exceeded the maximum allowable 3' matches, we can return, saying this primer length is invalid
            #print('primer size %d invalid, has %d 3' matches and score %f' % (len(seq_two), threeprimeMatches, totalScore))
            valid = False
            return(0, valid)
        
    totalScore = get_score_from_indices(matchLocations)
            
    return(totalScore, valid)

def check_for_dimers(sequence, forwardData, reverseData):
    """Check primers for the potential to form dimers with each other."""
    forward = zip(forwardData['score'], forwardData['length']) # Form a list of tuples of primer size and score.
    reverse = zip(reverseData['score'], reverseData['length'])
    
    minScore = max(forward[0]) * max(reverse[0])
    reverseSeq = str(sequence)[::-1]
    
    for f in forward: # Iterate over every valid primer combination and select the pair with the least potential to form dimers.
        fPrimer = str(sequence.complement())[0:f[1] - 1]

        for r in reverse:
            rPrimer = reverseSeq[0:r[1] - 1]
            alignment = pairwise2.align.localms(fPrimer, rPrimer, 2, -1, -.5, -.1, one_alignment_only=False)
            
            for align1, align2, score, begin, end in alignment:
                matchLocations, baseLocations = get_indices_from_alignment(align1, align2, debug=False)
                score = get_score_from_indices(matchLocations, debug=False)
                
                if score < minScore:
                    minScore = score
                    bestPairScores = (f[0], r[0])
                    bestPairSizes = (f[1], r[1])
    
    return bestPairSizes, bestPairScores

def get_score_from_sequences(seqOne, seqTwo, returnMin=False, oneAlignment=False, debug=False):
    """Returns the score of bonding two given sequences. If more than one 
    alignment is present, returns the maximum score or minimum if returnMin is true."""
    alignment = pairwise2.align.localms(seqOne, seqTwo, 2, -1, -.5, -.1, one_alignment_only=oneAlignment)
    count = 0
    maxScore = 0
    minScore = 0
    
    for align1, align2, score, begin, end in alignment:
        matchLocations, baseLocations = get_indices_from_alignment(align1, align2, debug=False)
        score = get_score_from_indices(matchLocations, debug=False)
        
        if debug:
            print("\n" + align1 + "\n" + align2 + "\n")
        
        if count == 0 or score > maxScore:
            maxScore = score
            
        if count == 0 or score < minScore:
            minScore = score
    
    if returnMin:
        return minScore
    
    else:
        return maxScore

def get_score_from_indices(matchLocations, debug=False):
    """Given a list of indices where two sequences match, gives the score of their alignment."""
    totalScore = 0.0
    
    for i in range(len(matchLocations)):
        totalScore += 2 # for every location where there is a match between two sequences, we add 2 to the score
        if i < len(matchLocations) - 1: # subtract .5 from the score for every gap, and .1 for each additional base pair in the gap
            diff = matchLocations[i + 1] - matchLocations[i] # diff is the size of the gap between our current nucleotide match and the next one
            if diff > 1:
                totalScore -= (diff - 2) * .1 + .5
    
    if debug:
        print("totalScore = " + str(totalScore))
    
    return totalScore

def get_threeprime_matches(matchLocations, threeprimeStart, threeprimeLength, maxThreeprimeGap, debug=False):
    """Returns the number of matches in the 3' region of a primer given a list of matching locations."""
    
    mostThreeprimeMatches = 0
    threeprimeMatches = 0
    threeprimeStartIndex = False
    
    for j in range(threeprimeStart, max(matchLocations) + 1): # Handle the case where the first base of the 3' region is not a match.
        if matchLocations.count(threeprimeStart) != 0:
            threeprimeStartIndex = matchLocations.index(threeprimeStart)
            break
        
        else:
            threeprimeStart = j + 1
            
    if not(threeprimeStartIndex):
        return 0
    
    for i in range(threeprimeStartIndex, len(matchLocations)):
        threeprimeMatches += 1
        mostThreeprimeMatches = max(mostThreeprimeMatches, threeprimeMatches)
        
        if i < len(matchLocations) - 1:
            if matchLocations[i + 1] - matchLocations[i] > maxThreeprimeGap + 1:
                threeprimeMatches = 0
                
    if debug:
        print("threeprimeMatches: " + str(mostThreeprimeMatches))
            
    return mostThreeprimeMatches

def get_indices_from_alignment(align1, align2, debug=False):
    """Returns the indices of matching nucleotides in a given alignment."""
    matchLocations = ''
    
    
    basesCounted = 0
    threeprimeMatches = 0
    mostThreeprimeMatches = 0
    
    align1Matches = [match.start() for match in (re.finditer(r'\w', str(align1)))] # Use regexes to find the indices of all characters that are not dashes.
    align2Matches = [match.start() for match in (re.finditer(r'\w', str(align2)))]

    matchLocations = [val for val in align1Matches if val in align2Matches] # The intersection of the above two lists (indices of all matches).
                
    if debug:
        print('\n' + align1)
        print(align2 + '\n\n')
        
    return matchLocations, align1Matches

def parse_sequence(sequence, seq_type):
    """Return a Seq object given a filename or string and a sequence type."""
    if seq_type == 'fasta': # if user supplies .fasta file, parse it using SeqIO
        
        handle = open(sequence)
        
        record_iterator = SeqIO.parse(handle, 'fasta')        
        rec_one = record_iterator.next()
    
        handle.close()
        
        seq_one = (rec_one.seq).upper() # sequence one is the target sequence
        
    elif seq_type == 'string': # if the user manually inputs a sequence of characters
        seq_one = Seq(sequence.upper(), generic_dna)
        
    else:
        #print('file type not recognized')
        seq_one = ''
        #leave for other types of sequences
        
    return seq_one

def get_sub_sequences(numSteps, seq_one, seq_two):
    """Split the template strand into a given number of subsequences and return a list of strings."""
    subSeqs = []
    
    if numSteps == 0 or len(seq_one) / 2 < len(seq_two):
        return [seq_one]
    
    else:
        subSeqs.append(seq_one)
        subSeqs.extend(get_sub_sequences(numSteps - 1, seq_one[0:len(seq_one) / 2], seq_two))
        subSeqs.extend(get_sub_sequences(numSteps - 1, seq_one[len(seq_one) / 2:], seq_two))
        return subSeqs
