------------------------------------------------------------------------
		SMETANA network alignment algorithm
------------------------------------------------------------------------
SMETANA finds the multiple alignment of a set of (biological) networks

SMETANA is written in MATLAB.

Release Note:
    10.31.2018
	1. A bug fixed in input networks read code
	2. Optimizing network read code  

Matlab Implementaion:
 Usage:
	alignment=SMETANA(Net_id_list,input_folder,id_flag,out_file)

 
   input   Net_id_list     the input networks ids list
           input_folder    the foder, where input files are placed
           id_flag         This flag is used to expedite the reading
                           process of inout files. The flag that indicates
                           whether the nodes are in numeric format
                           (id_flag=1) or not (id_flag=0). If nodes are in
                           the formatof "species id+number", the id_flag
                           should be 1, unless it is zero. For instance if
                           nodes are as a0,a1,b1,b3,... the id_flag is 1.
           out_file        output file address where the alignment result
                           will be written there
    output  alignment       Alignment result including the obtained aligned
                           nodes

Example:
        alignment=SMETANA({'dm','hs','sc'},'test',1,'output.txt')
        alignment=SMETANA({'a','b','c'},'test',1,'output.txt')

   Input files format:
       The Net_id_list provides the list of networks (spicies) names. For example
       here we have three species named as 'a', 'b', and 'c'. For these
       specises we need the following files:

       - Network files: a.net, b.net, c.net
               These are tab-separated files that list the (undirected)
               interactions in each network. For instanse, a.net may be as
               follows:

               a1	a2
               a3	a1
               a4	a2
               a2	a3

               Network files can also include the edge weights. In tha
               case we may have a third column:

               a1	a2  0.5
               a3	a1  0.2
               a4	a2  0.8
               a2	a3  0.9

       - Similarity scores files: a-b.sim, a-c.sim, b-c.sim

               These files consist of the list of similarity scores
               between nodes of different species. We have used BLAST
               bit-score for our tests, but other scores such as
               log(E-value) can also be used.

               IMPORTANT note 1:
               The similarity filename should have the species names in
               lexicographic order, i.e., a-b.sim is expected, not
               b-a.sim.

               IMPORTANT note 2:
               For a-b.sim file, the first column should contain node from
               speices a and the second column should contain IDs from
               species b.

               As an example, a-b.sim may be as follows:

               a1	b1	153
               a1	b3	55
               a1	b7	49
               a2	b3	444
               a3	b3	211
               a3	b4	122
               a4	b5	251
               a4	b8	71


    Output files format:
       The Output files includes the list of aligned nodes. For instance:

               a4 b5 c4
               a1 b1 b7 c1
               a2 b3 c2
               a3 b4 c3

 For more information on the algorithms, please see:

 Sahraeian, S.M.E., Yoon, B.J. (2013), SMETANA: Accurate and Scalable
 Algorithm for Probabilistic Alignment of Large-Scale Biological Networks,
 submitted.

 By Sayed Mohammad Ebrahim Sahraeian and Byung-Jun Yoon
 Contact: bjyoon@ece.tamu.edu
