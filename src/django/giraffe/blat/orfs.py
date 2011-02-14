
import models
from Bio.Seq import Seq

trans_table = 1 # standard translation table
min_protein_len = 200

def detect_orfs(sequence_object):
    features = []
    # each should be Sequence_Feature_Annotated()

    seq = Seq(sequence_object.sequence)
    seq_len = len(seq)

    for strand,nuc in [(+1,seq), (-1,seq.reverse_complement())]:
        for frame in range(3):
            #print 'strand '+str(strand)+' frame '+str(frame)

            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            # go through the translation one by one, so we know where
            # the ORFs start and end
            while aa_start < trans_len:
                #print 'start '+str(aa_start)
                aa_end = trans.find("*", aa_start)
                start_codon = trans.find('M', aa_start)
                has_stop = 1
                
                if aa_end == -1:
                    aa_end = trans_len
                    has_stop = 0

                # is there a start codon?
                if start_codon == -1 or start_codon > aa_end:
                    assert(aa_end != -1)
                    aa_start = aa_end+1
                    continue
                
                if aa_end-start_codon >= min_protein_len:
                    #print 'found '+trans[aa_start:aa_end]

                    # the following start and end need to start with
                    # 1, not 0.
                    if strand == 1:
                        start = frame+start_codon*3+1
                        end = min(seq_len,frame+aa_end*3+has_stop*3)
                    else:
                        start = seq_len-frame-aa_end*3-has_stop*3+1
                        end = seq_len-frame-start_codon*3

                    f = models.Sequence_Feature_Annotated()
                    f.sequence = sequence_object
                    f.feature_name = 'ORF frame '+str(frame+1)
                    f.feature_type_id = models.Feature_Type.ORF
                    f.orf_frame = frame
                    f.start = start
                    f.end = end
                    if strand == 1:
                        f.clockwise = True
                    else:
                        f.clockwise = False
                    features.append(f)
                    print str(f.to_dict())

                aa_start = aa_end+1

    sequence_object.clear_annotated_features(models.Feature_Type.ORF)
    for feature in features:
        feature.save()

