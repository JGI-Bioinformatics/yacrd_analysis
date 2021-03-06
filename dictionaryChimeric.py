from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
chimericExciseAndSplit = True
notCoveredExciseAndSplit = False
# read through yacrd results and store them in a diction; key is the read name, value is the rest of the information for that read (i.e. type such as 'Chimeric', length of read, and the zero coverage region information)
#yacrdFile = open('smallYacrd.txt')
yacrdFile = open('/global/projectb/scratch/ashleigh/chimeric/X0167/sequel/reads.yacrd')
yacrdEntries = {}
nameIndex = 1
for line in yacrdFile:
	line = line.strip().split()
	value = [line[0], line[2], line[3]]
	yacrdEntries[line[nameIndex]] = value
# read through fastq file, look up each read name in yacrdEntries to see if it is Chimeric or another type that needs changes
recordsToWrite = []
for record in SeqIO.parse('/global/projectb/scratch/ashleigh/chimeric/X0167/sequel/sequel.fastq', 'fastq'):
	identifier = record.id
	if identifier in yacrdEntries:
		regions = ''
		if ';' in yacrdEntries[identifier][2]:
			regions = yacrdEntries[identifier][2].split(';')
		else:
			reg = yacrdEntries[identifier][2]
			regions = []
			regions.append(reg)
		# add check here to see if length of not covered valid region is more than 25% of read, if so then go into below if
		if yacrdEntries[identifier][0] == 'Chimeric' or  yacrdEntries[identifier][0] == 'Not_covered':
			#print(record.name)
			numSubset = 0
			newReads = []
			if chimericExciseAndSplit:
				startIndex = 0
				indices = []
				allSegmentInfo = []
				#print(len(regions))
				for segment in regions:
					segmentInfo = segment.split(',')
					if len(segmentInfo) != 3:
						break
					#print(segmentInfo)
					allSegmentInfo.append(segmentInfo)
				if len(allSegmentInfo) == 0:
					continue
				if int(allSegmentInfo[0][1]) != 0:
					#means first chimeric region doesn't start at the beginning
					l = [0, int(allSegmentInfo[0][1]) - 1]
					indices.append(l)
				for i in range(len(allSegmentInfo)-1):
					l = [int(allSegmentInfo[i][2])+1, int(allSegmentInfo[i+1][1])-1]
					indices.append(l)
				index = len(allSegmentInfo)-1
				if int(allSegmentInfo[index][2]) != len(record.seq):
					l = [int(allSegmentInfo[index][2]) + 1, len(record.seq)]
					indices.append(l)
				if yacrdEntries[identifier][0] == 'Not_covered':
					sectionsLengthSum = 0
					for section in indices:
						sectionsLengthSum = sectionsLengthSum + (int(section[1]) - int(section[0]))
					if sectionsLengthSum < (.25 * int(yacrdEntries[identifier][1])):
						#print('entry ', identifier, ' does not meet 25% minimum coverage so it will not be added to the new fastq file.')
						continue
				for section in indices:
					newRead = record.seq[int(section[0]):int(section[1])+1]
					newQuality = record.letter_annotations['phred_quality'][int(section[0]):int(section[1])+1]
					name = identifier + '_' + str(numSubset)
					newRecord = SeqRecord(Seq(str(newRead), record.seq.alphabet), id = record.id + '_' + str(numSubset), description = record.description + '_' + str(numSubset), name = name)
					newRecord.letter_annotations['phred_quality'] = newQuality
					recordsToWrite.append(newRecord)
					#print(name)
					numSubset = numSubset + 1
	else :
		# write the original read directly to the file as is
		recordsToWrite.append(record)
SeqIO.write(recordsToWrite, '/global/projectb/scratch/ashleigh/chimeric/X0167/sequel/dictionaryOutput.fastq', 'fastq') #this one is for doing nothing with 'Not_covered'
