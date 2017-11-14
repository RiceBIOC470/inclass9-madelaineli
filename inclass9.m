% Inclass assignment 9
% GB comments
1) 100
2) 100
3) 100 Ratio in matches is miscalculated but will give credit because the question may be stated awkward. 
4) 100
Overall: 100

% The accession number for human NOTCH1 mRNA is AF308602
% 1. Read the information from this entry into matlab
accession = 'AF308602';
gb_data = getgenbank(accession);

% 2. Write code that runs a blast query on the first 500 base pairs of this
% gene against the refseq_rna database
seq_begin = gb_data.Sequence(1:500);
[requestID, requestTime] = blastncbi(seq_begin,'blastn','Database','refseq_rna');
blast_data = getblast(requestID,'WaitTime',requestTime*1.5);

% 3. Find the three highest scoring hits from other species and identify
% the length of the alignment and fraction of matches/mismatches. 

first = blast_data.Hits(2);
first_score = first.HSPs.Score
first_length = first.Length
first_align = first.HSPs.Alignment;
first_identity = first.HSPs.Identities;
first_match = first_identity.Match;
first_fraction = first_match/(500-first_match)

second = blast_data.Hits(6);
second_score = second.HSPs.Score
second_length = second.Length
second_align = second.HSPs.Alignment;
second_identity = second.HSPs.Identities;
second_match = second_identity.Match;
second_fraction = second_match/(500-second_match)

third = blast_data.Hits(7);
third_score = third.HSPs.Score
third_length = third.Length
third_align = third.HSPs.Alignment;
third_identity = third.HSPs.Identities;
third_match = third_identity.Match;
third_fraction = third_match/(500-third_match)

% 4. Run the same query against the database est_human. Comment on the
% sequences that you find. 

[requestID_2, requestTime_2] = blastncbi(seq_begin,'blastn','Database','est_human');
blast_data_2 = getblast(requestID_2,'WaitTime',requestTime_2*1.5);
highest = blast_data_2.Hits(1).HSPs;
second_high = blast_data_2.Hits(2).HSPs;
% When run against the database est_human, the length of alignment tend to 
% be much shorter than against refseq_rna. The score also tend to be
% smaller than against refseq_rna. For example, the highest scoring
% alignment against est_human has a length of 279, while the second
% hightest scoring alignment jumps straight to 35, which is very different
% from the query against refseq_rna. The highest scoring alignment in
% est_human database was against homo sapiens, (99% matching), while the second highest
% aligment was just very short sequence that just happen to align with the
% NOTCH1 given.
