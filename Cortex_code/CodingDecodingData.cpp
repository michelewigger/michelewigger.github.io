#include "CodingDecodingData.h"

namespace caching{

// A DFS based recursive function that returns true if a matching
// for vertex u is possible
bool bpm(vector<vector<bool> > bpGraph, int u, bool seen[], int matchR[])
{
    //Initialzation
    int N = bpGraph[0].size();
    // Try every job one by one
    for (int v = 0; v < N; v++)
    {
        // If applicant u is interested in job v and v is not visited
        if (bpGraph[u][v] && !seen[v])
        {
            // Mark v as visited
            seen[v] = true;

            // If job 'v' is not assigned to an applicant OR previously assigned
            // applicant for job v (which is matchR[v]) has an alternate job available.
            // Since v is marked as visited in the above line, matchR[v] in the following
            // recursive call will not get job 'v' again
            if (matchR[v] < 0 || bpm(bpGraph, matchR[v], seen, matchR))
            {
                matchR[v] = u;
                return true;
            }
        }
    }
    return false;
}

// Returns maximum number of matching from M to N
vector<vector<bool> > maxBPM(vector<vector<bool> > bpGraph)
{
    // An array to keep track of the applicants assigned to jobs.
    // The value of matchR[i] is the applicant number assigned to job i,
    // the value -1 indicates nobody is assigned.
    int M = bpGraph.size();
    int N = bpGraph[0].size();
    int matchR[N];

    // Initially all jobs are available
    for(int i=0; i<N; i++)
        matchR[i] = -1;

    // Count of jobs assigned to applicants
    int result = 0;
    for (int u = 0; u < M; u++)
    {
        // Mark all jobs as not seen for next applicant.
        bool seen[N];
        for(int i=0; i<N; i++)
            seen[i] = false;

        // Find if the applicant 'u' can get a job
        if (bpm(bpGraph, u, seen, matchR)){
            result++;
        }

    }
    /*cout << endl;
    for(int k=0; k<N; k++)
        cout << matchR[k] << " ";
    cout << endl;
    cout << "Nb jobs are: " << result << endl;*/

    //Graph to return based on matchR array
    vector<vector<bool>> r(M, (vector<bool> (N,false)));
    //int cc=0;
    for(int k=0; k<N; k++){
        if(matchR[k]  != -1){
            r[matchR[k]][k] = true;
            //cc++;
        }
    }

    //test debug
    /*for(int i=0; i<M; i++){
        for(int j=0; j < N; j++){
            cout << r[i][j] << " ";
        }
        cout << endl;
    }
    cout << "Total nb of jobs: " << cc << endl;*/

    return r;
}

// Fill in the Q vector of requested files by strong users
void gen_rand_request_zipf(int m_files, double alpha, vector<int> &Q){

    int idx_rand;
    int nb_users = Q.size();
    unsigned long int randSeed = 0;
    double probs_vec [m_files];
    double sum=0;
    /******************Files popularity generation*********************/
    for (int i=0; i<m_files; i++)
    {
        double j = i+1;
        probs_vec[i] = (1/pow(j,alpha));
        sum += probs_vec[i];
    }
    for(int i=0; i<m_files; i++)
        probs_vec[i] = probs_vec[i] / sum;

    /******************File index to request generation***************/
    gsl_rng *r;
    gsl_ran_discrete_t *rand_disc;

    srand(time(NULL));                    /* initialization for rand() */
    randSeed = rand();                    /* returns a non-negative integer */

    gsl_rng_env_setup();
    gsl_rng_default_seed = randSeed;

    r = gsl_rng_alloc(gsl_rng_default);

    /*For each user generate a random request*/
    for (int i=0; i<nb_users; i++){
        rand_disc = gsl_ran_discrete_preproc (m_files, probs_vec);
        idx_rand =  gsl_ran_discrete (r, rand_disc);

        if (idx_rand >= m_files || idx_rand < 0){
            printf("\nError: Random Index.\n");
            exit(0);
        }

        Q[i] = idx_rand;
    }
}
/*! This function creates the bipartite graph and apply the matching between nodes
    in order to code the strong adn weak packets together.
    It also reads the packets requested from the strong users
    @param  coloring is an array that contains the colors of all the nodes in conflict graph
    @param  n_col is defined as the total number of colors in the conflict graph that
            also represents the total number of transmitted packets
    @param  nodi is an array that contains all the nodes of the conflict graph
    @param  nb_nodes is the total number of nodes
    @param  nb_strg is the number of strong users
    @param  data is a structure thatcontains the configuration loaded from the environment
            it contains the nb of users, nb of files in te library, nb of chunks of each fle, cache matrix,
            vector of request, and packet to transmit
    @param  header_data is a 2d array that is given empty and is filled by the headers of the strong packets
    @param  G_edges is an adjacent matrix between the strong nodes and the weak nodes. It indicates the edges of the graph
    \return a 2d vector that contains the strong packets to be polarly coded,
            and fill in the header data of these packet, in addtion to the edges matrix of the graph
*/

vector<vector<char>> MaxBipartiteGraph(int *coloring, int n_col, nodo *nodi, int nb_nodes, int nb_strg, data_matrix data,
    header_transmission **header_data, std::vector<std::vector<bool> > &G_edges_vec){  //bool*** G_edges

    //Combine_nodes is a 2d vector, the rows are the number of transmission,
    //and at the columns are the number of users in each transmitted packets
    vector<vector<nodo>> node_w;
    vector<nodo> node_s;
    int nb_weak = data.n_utenti;
    int nb_chunks = data.b_chunks;
    int nb_user =  nb_weak + nb_strg; //nb of weak users
    int m_files = data.m_files;
    nodo strongNode;
    int userid, fileid, chunkid, size_strg, size_weak;
    double alpha = 0.0;
    //edges matrices
    vector<vector<bool>> Matrix_Adj;
    //int*** u_cache = data.Ind;
    bool edge;


    /*--------------------CORTEXLAB-------------------*/
    //Path to extract data chunks
    string pathFiles = "../repository"; //"/CachingFile/repository"
    vector<string> files; //all files of the reository

    int debug =1;


    //loop for the rows ; the number of transmission
    for (int i = 0; i < n_col; i++){
        int color = i + 1;
        vector<nodo> v_i;
        //loop for the users in each packet
        for (int j = 0; j < nb_nodes; j++){
            if (coloring[j] == color)
                v_i.push_back(nodi[j]);
        }
        node_w.push_back(v_i);
    }

    /*************File index to be requested generator for Strong users*******/
    vector<int> Q(nb_strg, 0);
    for (int i=0; i<nb_strg; i++)
        gen_rand_request_zipf(m_files, alpha, Q);
    if(debug){
        cout << "Files requested by strong users: ";
        for (int i = 0; i < nb_strg; ++i)
            cout << Q[i] << ", ";
        cout << endl;
    }

    //Create the nodes corresponding to the strong users request
    int ct = 0;
    strongNode.id_chunck.push_back(0);
    for (int i = nb_weak; i < nb_user; i++) {
        for (int j = 0; j < nb_chunks; j++) {
            ct ++;
            strongNode.id_utente = i;
            strongNode.id_file   = Q[i-nb_weak];
            strongNode.id_chunck[0] = j;
            strongNode.id        = ct;
            strongNode.degree    = 0;
            node_s.push_back(strongNode);
        }
    }
    size_strg = node_s.size();
    size_weak = node_w.size();

    /*if(debug){
        cout << "Strong node size = " << size_strg << endl;
        cout << "Weak node size = "   << size_weak << endl;
        for (int i = 0; i < size_strg; i++)
            cout << "id: " << node_s[i].id << ", idUser: " << node_s[i].id_utente
                << ", idFile: " << node_s[i].id_file << ", idChunk: " << node_s[i].id_chunck << endl;

    }*/

    /**********Build the Bipartite graph****************/

    //Initialize the Bpartite graph adjacent matrix
    /*Matrix_Adj = (bool**) malloc(node_s.size()*sizeof(bool* ));
    for(unsigned int i=0; i<node_s.size(); i++)
        Matrix_Adj[i] = (bool*) malloc(node_w.size()*sizeof(bool));*/
    for (int i = 0; i < size_weak; i++){
        vector<bool> temp(size_strg, false);
        /*for (int j = 0; j < size_strg; j++)
            temp.push_back(false);*/
        Matrix_Adj.push_back(temp);
    }

    if(debug){
        cout << "adjacent matrix size: " << Matrix_Adj.size() << endl;
    }

    //Check for every strong packet, the weak nodes that can be combined with
    for (int i = 0; i < size_strg; i++) {
        chunkid = node_s[i].id_chunck[0];
        fileid = node_s[i].id_file;
        ct=0; //counter for possible combination of strong and weak
        /*if(debug)
            cout << "File ID, Chunck ID: " << fileid << ", " << chunkid << endl;*/

        //Go over all the weak nodes for one strong node
        for(int j=0; j<size_weak; j++){
            edge = true;
            // for one weak packet and one strong node, check if the requested packet
            //of the strong user is included in the cache of all weak users
            /*if (i==0){
                cout << ct << " mixed nodes = " << node_w[j].size();
                //ct++;
            }*/
            //IN CASE NO XOR APPLIED, THE  USERID WILL BE THE SAME FOR ALL u
            for(unsigned int u=0; u < node_w[j].size(); u++){
                userid = node_w[j][u].id_utente;

                if(data.Ind[userid][fileid][chunkid] == 0)
                    edge = false;
            }
            /*if (i==0)
                cout << endl;*/
            if(edge){
                Matrix_Adj[j][i] = true;
                ct++;
            }
            else
                Matrix_Adj[j][i] = false;
        }

    }

    if(debug){
        int cc = 0;
        for (int i = 0; i < size_weak; i++) {
            for(int j=0; j< size_strg; j++){
                if(Matrix_Adj[i][j]==true){
                    cc++;
                }
            }
        }
        cout << "Nb of edges = " << cc << endl;
        //cout << "portion of combined packet, 1 out of " << (size_weak*size_strg/cc) << endl;
    }


    //Apply te maximum bipartite matching algorithm
    /*vector<vector<bool>> G_edges_vec;
    bool** Edges_tmp;*/
    int sum=0;
    G_edges_vec = maxBPM(Matrix_Adj);
    for(int i=0; i<G_edges_vec.size(); i++){
        for (int j = 0; j < G_edges_vec[i].size(); ++j)
        {
            if(G_edges_vec[i][j]){
                sum++;
                //cout << j << ", ";
            }
        }
    }
    cout << "Number of Edges after Max Bipartite Algo: " << sum << endl;

    /*
    //Allocation of the pointer G-edges that will be modified by reference
    Edges_tmp =  (bool**) malloc(size_weak * sizeof(bool*));
    //check_memory_allocation_2D(Edges_tmp, "Allocation 2D Edges_tmp.");

    for(int i=0; i<size_weak; i++){
        Edges_tmp[i] = (bool*) malloc(size_strg * sizeof(bool));
        //check_memory_allocation_1D(Edges_tmp[i], "Allocation 1D Edges_tmp[i]");
    }
    //Populate the pointer from the reduced matrix G_edges_vec
    for(int i=0; i<Matrix_Adj.size(); i++){
        for(int j=0; j<Matrix_Adj[i].size(); j++){
            Edges_tmp[i][j] = G_edges_vec[i][j];
        }
    }
    (* G_edges) = Edges_tmp;
    */
    /*if(debug){
        int cc=0; // This variable counts the total number of existing edges
        int cccc=0; // This variable counts the total nb of strong/weak packages that can be combined
        for(int i=0; i<size_weak; i++){
            int ccc=0; // counts te nb of edges for every weak packet
            for(int j=0; j<size_strg; j++){
                if(Edges_tmp[i][j]){
                    cc++;
                    ccc++;
                }
            }
            //cout << ccc << " ";
            if(ccc >= 1)
                cccc++; //increment by one if there are one or more edges
        }
        cout << "\nThe Nb of edges after filtering: " << cc << endl; //<< ". Possible combination = " << cccc << endl;
    }*/


    //Create the header of each packet
    header_transmission *header = new header_transmission[size_strg];
    if (!header){
        cout << endl << "MEMORY EXCEPTION: header memory allocation." << endl << endl;
        exit(0);
    }

    //Load the names of the files into an array
    files = getDirectoryFiles(pathFiles);
    //for(int i=0; i< files.size(); i++)
        //cout << files[i] << ", ";
    //cout << endl;

    if (files.empty()){
        cout << endl << "Repository folder is empty or it not exist." << endl << endl;
        exit(0);
    }

    /*Sort files by name*/
    sort(files.begin(), files.end());

    /*Find max size on each package*/
    double max_size = -INF;
    for (string f : files){
        int size_file = getFileSize(f);
        //cout << size_file << ", ";
        if (size_file > max_size){
            max_size = size_file;
        }
    }
    //cout << endl;
    int max_size_package = ceil(max_size / (double) nb_chunks);

    cout << "The max size package is: " << max_size << "/" << nb_chunks << " = " << max_size_package << endl;
    //cout << "The file size: " << max_size << endl;

    vector< vector<char> > strgusr_data(size_strg, vector<char>(max_size_package, 0));

    /*For each node create a chunk ready for transmission*/
    for (int i = 0; i < size_strg; i++){

        nodo strg_node = node_s[i];

        userid = strg_node.id_utente;
        fileid = strg_node.id_file;
        chunkid = strg_node.id_chunck[0];

        header[i].id_utenti.push_back(userid);
        header[i].id_files.push_back(fileid);
        header[i].id_chunks.push_back(chunkid);
        header[i].nb_chunks.push_back(1);

        /*Reading file */
        string file = files.at(fileid);
        ifstream is (file, ifstream::binary);
        //If the file is successfully opened
        if (is) {
            double size_file = getFileSize(file);
            int size_package = ceil(size_file / (double) nb_chunks);

            int begin_package = (chunkid * size_package);
            //cout << begin_package << ", ";

            if (begin_package > size_file){
                begin_package = size_file;
                size_package = 0;
            }else{
                if (begin_package + size_package > size_file){
                    size_package = abs(size_file - begin_package);
                }
            }

            if (size_package > 0){
                is.seekg(begin_package);

                char *buffer = new char[size_package];

                is.read(buffer, size_package);

                is.close();

                for (int k = 0; k < size_package; k++)
                    strgusr_data.at(i).at(k) = buffer[k];
                /*for(int k=0; k<20; k++)
                    cout << (int)strgusr_data.at(i).at(k) << ", ";
                cout << endl;*/
                //cout << size_package << ", ";
                delete[] buffer;
            }

            header[i].size_package.push_back(size_package);
        }else{
            cout << endl << "Error opening file in storage: " << file << endl;
            exit(0);
        }
    }

    (*header_data) = header;
    return strgusr_data;

}


/*! This function creates the coded data to be sent after being polarly coded strong and weak
    based on the adjacent matrix of the bipartte graph
    \param  weak_data is a matrix that contains in its rows the weak packages and in its columns the data of each package.
            The size is nb_weak_packages x size_package
    \param  strg_data is a matrix that contains in its rows the strong packages and in its columns the data of each package.
            The size is nb_strong_packages x size_package
    \param  G_edges is an adjacent matrix of the Bipartite Grapgh
    \param  hdr_weak is an array that contains the headers of the weak packets - 2d vector
    \param  hdr_strg ia an array that contains the headers of the strong packets - 2d vector.
            This can be reducd to 1d vector, however this would complicate the header combination (to be verified later)
*/
vector<vector<char>> codingDataPolar(vector<vector<char>> weak_data, vector<vector<char>> strg_data, vector<vector<int> > &data_bits,
    vector<vector<bool>> G_edges, header_transmission *hdr_weak, header_transmission *hdr_strg, vector<header_polar> &hX, const int N){

    bool DEBUG = false;
    int nodeweak_size = weak_data.size();
    int nodestrg_size = strg_data.size();
    int ct=0;

    /******************************** Final Packet Before Encoding *********************/
    vector<vector<char> > PC_data;
    //A vector that contains the IDs of strong packets successfully combined with the weak packets
    vector<int> comb_strg_pck;
    int j, cc, cl, cn, cw;
    //FOR TEST ONLY
    for (int i=0; i< weak_data.size(); i++)
        cout << weak_data.at(i).size() << ", ";
    cout << endl;
    for (int i=0; i< strg_data.size(); i++)
        cout << strg_data.at(i).size() << ", ";
    cout << endl;
    cout << "--- Coding Data Polarly ---" << endl;
    //For each weak packet, check if it can be combned with a strong one
    for(int i=0; i<nodeweak_size; i++){
        vector<char> vd;
        //vd = strg_data.at(i);
        //For each strong packet, check if the weak packet can be combined with
        for(j=0; j<nodestrg_size; j++){
            //If there is an edge between one strong and one weak, they are sent in one packet
            if(G_edges[i][j]){
                unsigned int len_w = weak_data.at(i).size();
                for(unsigned int k=0; k<len_w; k++){
                    vd.push_back(weak_data.at(i).at(k));
                    //if(i<5)
                    //    cout << (int) weak_data.at(i).at(k) << ", ";
                }
                //if(i<5)
                //    cout << endl;
                unsigned int len_s = strg_data.at(j).size();
                for(unsigned int k=0; k < len_s; k++){
                    vd.push_back(strg_data.at(j).at(k));
                }
                //cout << "[" << len_w << "," << len_s << "]";
                comb_strg_pck.push_back(j);
                //Construct the header
                header_polar hd; //tmp header
                cc = 0;
                for(unsigned int k=0; k<hdr_weak[i].id_utenti.size(); k++){
                    hd.id_utenti.push_back(hdr_weak[i].id_utenti[k]);
                    hd.id_files.push_back(hdr_weak[i].id_files[k]);
                    hd.nb_chunks.push_back(hdr_weak[i].nb_chunks[k]);
                    for(unsigned int l=0; l < hdr_weak[i].nb_chunks[k]; l++)
                        hd.id_chunks.push_back(hdr_weak[i].id_chunks[cc+l]);
                    cc+=hdr_weak[i].nb_chunks[k];
                }
                //cout << cc << ", ";

                hd.id_utenti.push_back(hdr_strg[j].id_utenti[0]);
                hd.id_files.push_back(hdr_strg[j].id_files[0]);
                hd.id_chunks.push_back(hdr_strg[j].id_chunks[0]);
                hd.nb_chunks.push_back(1);
                //cout << hdr_strg[j].id_chunks[0] << ", ";
                //hd.size_package.push_back(hdr_strg[j].size_package[0]);
                hd.strong = true;
                hd.weak = true;
                hX.push_back(hd);
                //There is no need to check the other strong packet since only one edge exists, then break
                break;
            }
        }
        //If strong and weak can be polarly coded together
        if(j<nodestrg_size){
            PC_data.push_back(vd);
        }else{ // If the weak cannot be polarly coded with another strong, then weak padded with zeros
            int len_w = weak_data.at(i).size();// SHOULD BE EQUAL TO THE LENGTH OF STRONG PACKET
            for(int k=0; k < len_w; k++){
                vd.push_back(weak_data.at(i).at(k));
            }

            int len_s = strg_data.at(0).size();
            for(int j=0; j<len_s;j++)
                vd.push_back(0);

            PC_data.push_back(vd);
            //Construct the header
            header_polar hd; //tmp header
            cc = 0;
            for(unsigned int k=0; k<hdr_weak[i].id_utenti.size(); k++){
                hd.id_utenti.push_back(hdr_weak[i].id_utenti[k]);
                hd.id_files.push_back(hdr_weak[i].id_files[k]);
                hd.nb_chunks.push_back(hdr_weak[i].nb_chunks[k]);
                for(unsigned int l=0; l < hdr_weak[i].nb_chunks[k]; l++)
                    hd.id_chunks.push_back(hdr_weak[i].id_chunks[cc+l]);
                cc+=hdr_weak[i].nb_chunks[k];

                //hd.size_package.push_back(hdr_weak[i].size_package[k]);
                //cout << hdr_weak[i].id_utenti[k] << ", ";
            }
            hd.strong = false;
            hd.weak = true;
            hX.push_back(hd);
        }

    }
    //Print out the nb packets and the combined ones
    cout << "Number of requested Strong Packets: " << nodestrg_size << endl;
    cout << "\nNumber of Packets - weak and weak with strong: " << PC_data.size() << endl;
    int NbStrgCombPack = comb_strg_pck.size();
    int RemainStrgPack = nodestrg_size - NbStrgCombPack;
    cout << "Number of Strong Packets combined with weak: " << NbStrgCombPack << endl;
    cout << "Minimum Rate Gain: " << (PC_data.size()+nodestrg_size)/((float) PC_data.size()+NbStrgCombPack) << endl;
    /*for(int k=0; k<NbStrgCombPack; k++)
        cout << comb_strg_pck[k] << " ";
    cout << endl;*/

    //Add all remaning strong packets at the end of the data
    //vector<vector<char>> vd;
    cc=0, cl=0, cn=0;
    /*for(int k=0; k<strg_data.size(); k++)
        cout << strg_data[k].size() << " ";
    cout << endl;*/

    //HERE WE ASSUME THAT THE LENGTH OF THE PACKETS IS FIX FOR ALL THE CHUNKS OVER ALL FILES
    //Construct the header -- tmp header
    header_polar hd;
    int vd_size = strg_data[0].size();
    vector<char> tmp (4*vd_size, 0);
    //cout << nodestrg_size << "\t";
    for(int i=0; i<nodestrg_size; i++){

        //int vd_size = strg_data[i].size();
        if(find(comb_strg_pck.begin(), comb_strg_pck.end(), i) == comb_strg_pck.end()){
            //cout << ct++ << " ";
            for(int j=0; j<vd_size;j++){
                tmp[cn] = strg_data.at(i).at(j);
                cn++;
                //tmp.push_back(strg_data.at(i).at(j));
            }

            //vd.push_back(tmp);
            cc++;
            //cout << hdr_strg[i].id_chunks[0] << ", ";
            cl++;//Count the packets non combined
            //header addition part
            hd.id_utenti.push_back(hdr_strg[i].id_utenti[0]);
            hd.id_files.push_back(hdr_strg[i].id_files[0]);
            hd.id_chunks.push_back(hdr_strg[i].id_chunks[0]);
            hd.nb_chunks.push_back(1);

            //hd.size_package.push_back(hdr_strg[i].size_package[0]);
            hd.strong = true;
            hd.weak = false;
        }
        if(cc == 4 || cl==RemainStrgPack){
            //Add counter to debug
            ct++;
            cn=0;
           
            PC_data.push_back(tmp);
            hX.push_back(hd);

            /********FOR TEST ONLY***********/
            if(DEBUG){
                ofstream data_out;
                data_out.open("testStrongUsers.txt",ios::app);
                for (int j = 0; j < tmp.size(); ++j)
                    data_out << tmp[j];
                data_out << endl << endl ;
                data_out.close();
            }

            //still needs to add the frozen bits
            cc = 0;
            tmp = vector<char> (4*vd_size,0);
            //re-initialize the tmp header
            hd.id_utenti = vector<unsigned int>();
            hd.id_files  = vector<unsigned int>();
            hd.id_chunks = vector<unsigned int>();
            //hd.size_package = vector<unsigned int>();
        }
        if(cl==RemainStrgPack)
            break;

    }
    cout << "The number of strong packet combined together is " << ct << endl;

    //TEST
    /*for(int mm=0; mm < hX.size(); mm++){
        //if(hX[mm].id_utenti.size()>1){
            for (unsigned int ll = 0; ll < hX[mm].id_utenti.size(); ll++)
                cout << "[" << ll << "," << hX[mm].id_utenti.at(ll) << "," << hX[mm].id_files.at(ll) << "]";
            cout << endl;
            for(int kk=0; kk < hX[mm].id_chunks.size(); kk++)
                cout << hX[mm].id_chunks[kk] << ", ";
            cout << endl;
        //}
    }*/
    

    /******************************** Polar Encoding *********************/
    //const int N = 2048;
    //Polar coding -- adding zeros frozen bits included
    std::vector<unsigned int> bb (8,0); //bits vector of each converted char

    //data in bits
    //data_bits is the vector the data to polar code but in bits and padded with zeros;
    //std::vector<std::vector<int> > data_bits (PC_data.size(), std::vector<int> (N,0));
    data_bits = vector<vector<int> > (PC_data.size(), vector<int> (N,0));
    for(unsigned int i=0; i<PC_data.size(); i++){
        unsigned int packetSize = PC_data[i].size();
        //cout << packetSize << ", ";
        for(unsigned int j=0; j< packetSize; j++){
            //if(i<5)
            //    cout << (int)PC_data[i][j] <<", ";
            bb = conv_char_to_bitsInt(PC_data[i][j]);
            for(int k=0; k<8; k++)
                data_bits[i][j*8+k] = bb[k];
        }
        //if(i<5)
        //    cout << endl << endl;

    }

    //cout << endl << data_bits.size() << ", " << N << endl;

    /*for(int i=0;i<data_bits.size(); i++){
        int index = index_find(hX[i].id_utenti, 2);
        if(index!=-1){
            int sum=0;
            for(int j=0; j< data_bits[i].size(); j++){
                //if(i<4)
                //    cout << data_bits[i][j];
                sum+=data_bits[i][j];
            }
            cout <<  sum << ", "; //"\nsum_bits:" <<
        }
    }*/

    //----------------Freeing vectors-----------------
    comb_strg_pck.clear();
    for (int i = 0; i < PC_data.size(); ++i)
        for (int j = 0; j < PC_data[i].size(); ++j)
            PC_data[i].clear();
    //-------------------------------------------------


    //Encode using polar encoder
    PC PC_w, PC_s;
    //bits_coded are the bits after being polarly coded
    std::vector<std::vector<int> > bits_coded (data_bits.size(), std::vector<int> (N,0));
    //bits_coded = vector<vector<int> > (data_bits.size(), std::vector<int> (N,0));
    const int size_chunk = 8*strg_data.at(0).size();
    int K_w = 8*weak_data.at(0).size();
    //int K_s = K_w + size_chunk;
    int K_s = size_chunk*3 + size_chunk;

    cout << "Coding parameters: " << K_w << ", " << K_s << ", " << size_chunk << endl;


    //Coding intialization
    double designSNRdb = 0;
    PC_w = initialize_PC(N,K_w);
    PC_s.initPC(N, K_s, designSNRdb);
    PC_s.setGenMatrix(PC_w.genMatrix);
    PC_s.setRn(PC_w.Rn);
    PC_s.setArrangedBits(PC_w.arragedBits);

    int info_w [K_w], frozen_w[data_bits.size()][N-K_w], frozen_s[N-K_s];

    // Same for both users
    int sentMessage[N], sentCodeword[N], sentSymbol[N];//sentSymbol[data_bits.size()][N];

    //Polar Encoding
    //cout << "sum of the bits after each encoding: " << endl;
    for (unsigned int k = 0; k < data_bits.size(); ++k)
    {
        for (int i=0; i<N; i++){
            //sentMessage[k][i] = rand()%2;
            sentMessage[i] = data_bits[k][i];
        }

        for (int i=0; i<K_w; i++)
            info_w[i] = sentMessage[i];
        for (int i=K_w; i<N; i++)
            frozen_w[k][i-K_w] = sentMessage[i];


        PC_w.encode(info_w, frozen_w[k], sentMessage, sentCodeword, sentSymbol);
        for (int i=0; i<N; i++){
            bits_coded[k][i] = (sentSymbol[i]>0) ? 1:0;
            //sentSymbol[k][i] = 2*bits_coded[k][i]-1;

        }
        /*if(k < 20)
        {
            int sum=0;

            for (int i = 0; i < N; ++i){
                //cout << bits_coded[k][i];
                sum = sum+bits_coded[k][i];
            }
            cout << sum << ", "; //"\nsum bits coded: " <<
            //cout << endl << endl;
        }*/

    }
    //cout << endl;
    /******************************** Polar Decoding - FOR TEST *********************/
    //Function to implement/use at the receiver -- IT IS HERE FOR TEST

    //Test decoding -- Verify that procedure until settlement of frozen bits in all cases is well implemented
    //Note that here we assume that data_bits is known where it is not in reality
    //In general, the receiver should first read the header to know if it is concerned by the packet
    //If it is concerned:
    //   weak receiver checks if its packet is combined with strong packet
    //   if yes,
    //      the receiver check the ID of the strong packet and read it from the cache and complete frozen bits
    //      then proceed to polar decoding
    //    if no,
    //      set frozen bits to zeros and proceed to polar decoding
    //   Now, the weak packet is polarly decoded, it proceeds to classical decoding (implemented)
    //   Strong receiver polarly decodes and extract strong packet (implemented)



    double llr_w[N], llr_s[N];
    float recSymbol_w[N], recSymbol_s[N];
    int recMessage_w[N], recCodeword_w[N], recMessage_s[N], recCodeword_s[N];

    int user_id = 4;    
    int coderate[5] = {1, 2, 2, 3, 4};
    K_w = coderate[user_id]*size_chunk;
    int frozen_bits[N-K_w];
    cout << K_w << endl;
    PC PC_w1;
    PC_w1.initPC(N, K_w, designSNRdb);
    PC_w1.setGenMatrix(PC_w.genMatrix);
    PC_w1.setRn(PC_w.Rn);
    PC_w1.setArrangedBits(PC_w.arragedBits);

    //============QPSK mapping=============================
    cout << "------QPSK Mapping ------" << endl;;
    int Symb_size = N/2;
    vector<vector<gr_complex> > tx_Symb;
    gr_complex sentSymbol_qpsk[bits_coded.size()][Symb_size], recSymbol_qpsk[Symb_size];

    tx_Symb = BitsToQPSKSymb(bits_coded);

    for (int i = 0; i < tx_Symb.size(); ++i){
        for (int j = 0; j < tx_Symb[i].size(); ++j){
            sentSymbol_qpsk[i][j] = tx_Symb[i][j];
        }
    }
    cout << "------QPSK Mapping OK ------" << endl;
    
    //=====================================================
    ofstream data_out;
    float snrmin = 0, snrmax = 6, SNR_u;
    vector<double> ber_total;
    int index, NbChnl=200;
    int cwErrorSum_w, bitErrorSum_w, error_w;

    //Channel parameters
    const int N_carr = 64, log2N = 6;
    int Nb_packets = tx_Symb.size();
    int t[] = {0, 2, 4};
    double pow_delay[] = {0.5, 0.3, 0.2};
    int OFDM_Symb = ceil(N/(2*N_carr));
    int L=3, Ls=8, fd=40, T=80*pow(10,-6); 
    gr_complex **a;
    gr_complex h[N_carr],f[N_carr],yf[N_carr],in_ch[N_carr],out_ch[N_carr],out_ch_p[N_carr];
    gr_complex **hp;
    gr_complex noise, nh;
    //double snr_m=0;

    hp = (gr_complex **) malloc (bits_coded.size() * sizeof(gr_complex *));
    check_memory_complex_allocation_2D(hp, "Allocation 2D Probs.");

    for (unsigned int i=0; i<bits_coded.size(); i++){
        hp[i] = (gr_complex *) malloc (Symb_size * sizeof(gr_complex));
        check_memory_complex_allocation_1D(hp[i], "Allocation 1D Probs.");
    }

    data_out.open("Results.txt",ios::app);
    data_out << "\nuser: " << user_id << endl;
    
    /*default_random_engine gen( rand() );
    normal_distribution<float> n(0,1);
    normal_distribution<float> chn(0,10);
    double h_norm=0;
    for (int i = 0; i < 1000000; ++i){
        snr_m += pow(abs(divide_complex(gr_complex(n(gen),n(gen)),gr_complex(chn(gen),chn(gen)))),2);
    }
    snr_m = (double) 1000000/snr_m;*/

    //----------------Freeing vectors-----------------
    for (int i = 0; i < tx_Symb.size(); ++i)
        for (int j = 0; j < tx_Symb[i].size(); ++j)
            tx_Symb[i].clear();
    //-------------------------------------------------

    cout << "---------- Loop SNR --------" << endl;
    for(SNR_u = snrmin; SNR_u <=snrmax; SNR_u+=0.5)
    {
        float variance_u = pow(10,-(SNR_u/10));
        float sqrtVariance_u = sqrt(variance_u);
        //cout << "SNR: " << SNR_u << ", " << sqrtVariance_u*(var_vec[SNR_u-snrmin]) << endl ;
        cout << "SNR: " << SNR_u << ", " << sqrtVariance_u << endl ;

        //noise generator function
        //default_random_engine generator( rand() );
        //normal_distribution<double> distribution(0,sqrtVariance_u);

        double ber =0;
        for (int ch = 0; ch < NbChnl; ++ch)
        {
            
            cwErrorSum_w = 0; bitErrorSum_w = 0;
            //===================Channel generation=====================================
            /*a = generate_ch_param(pow_delay, OFDM_Symb*Nb_packets, N_carr, Ls, fd, T, L);

            for (int i = 0; i < Nb_packets; ++i) //
            {
                for (int j = 0; j < OFDM_Symb; ++j)
                {
                    // Forward transform
                    memcpy((void *) &f[0], (void *) &sentSymbol_qpsk[i][j*N_carr], sizeof(gr_complex)*N_carr);
                    //Invert FFT for OFDM multiplexing
                    iFFT( f, in_ch, log2N );          
                    //generate attenuation vector a1,a2,...,aL for one OFDM symb i.e. length N_carr
                    dyn_chnl(out_ch_p, in_ch, a[i*OFDM_Symb+j], N_carr, L, t);
                    //PC_w1.noise_gc(out_ch_p, out_ch, sqrtVariance_u,N_carr);//

                    //FFT of noiseless received signal for channel estimation
                    if(j==0){
                        FFT ( out_ch_p, yf, log2N );   //print( "\nTransform:", f, log2N );
                        //Channel estimation and equalization
                        float tot_norm =0;
                        for (int k = 0; k < N_carr; ++k){
                            h[k]= yf[k]/f[k];
                            tot_norm += pow(abs(h[k]),2);
                        }
                        tot_norm /= N_carr;
                        //cout << tot_norm << ",";
                        //for (int k = 0; k < N_carr; ++k){
                        //    h[k] /= sqrt(tot_norm);  
                        //}
                        
                    }
                    //FFT for OFDM demultiplexing
                    FFT ( out_ch, yf, log2N );   //print( "\nTransform:", f, log2N );
                    //Equalization
                    for (int k = 0; k < N_carr; ++k){
                        //yf[k]= yf[k]/h[k];
                        hp[i][j*N_carr+k] = h[k];
                    }

                    //for(int k=0; k<N_carr; ++k){
                        //noise = gr_complex(distribution(generator)/sqrt(2),distribution(generator)/sqrt(2));
                        //nh = divide_complex(noise, h[k]);
                        //yf[k] = yf[k];// + (nh);//
                        //rx_Symbol_qpsk[i][j*N_carr+k] = yf[k] + noise;
                        //cout << yf[k] << ", " << nh << ", ";
                    //}
                    //cout << endl << endl;

                    //REBUILD THE RECEIVED DATA SYMBOLS FROM yt
                    //memcpy((void *) &rx_Symbol_qpsk[i][j*N_carr], (void *) &yf[0], sizeof(gr_complex)*N_carr);
                            
                }
            }

            int len_a = OFDM_Symb*Nb_packets;
            for (int i = 0; i < len_a; ++i)
                delete [] a[i];*/
            //===================End Channel============================================

            for (int k = 0; k < bits_coded.size(); ++k){//bits_coded.size()
                error_w =0;

                //==========packet for which user===========
                cc=0;
                index = index_find(hX[k].id_utenti, user_id);
                if(index != -1){
                    //cout << k << ", " << K_w << endl;
                    
                    for (int i = 0; i < N-K_w; ++i){
                        frozen_bits[i] = data_bits[k][K_w+i];
                    }

                    //PC_w.noise(sentSymbol[k], recSymbol_w, sqrtVariance_u);
                    //PC_w.computeLLR(llr_w,recSymbol_w,variance_u);
                    PC_w1.noise_gc(sentSymbol_qpsk[k], recSymbol_qpsk, sqrtVariance_u,N/2);
                    //PC_w1.noise_hgc(sentSymbol_qpsk[k], recSymbol_qpsk, hp[k], sqrtVariance_u,N/2);
                    PC_w1.computeLLR_qpsk(llr_w,recSymbol_qpsk,variance_u);
                    PC_w1.SC(recMessage_w, recCodeword_w, llr_w, frozen_bits);
                    recMessage_w[0] = (recMessage_w[0])?0:1;
                    //Test for errors
                    for (int i=0; i<K_w; i++)
                        if (recMessage_w[i] != data_bits[k][i])
                            error_w = error_w + 1;

                    if (error_w != 0)
                        cwErrorSum_w = cwErrorSum_w + 1;

                    bitErrorSum_w = bitErrorSum_w + error_w;
                    cc++;
                }
            }
            cout << "Error per bit = " << bitErrorSum_w << endl;
            //cout << "SNR measured is " << snr_m << endl;
            ber += (double)bitErrorSum_w/(K_w*cc);
        }
        cout << "Average BER = " << (double)ber/NbChnl << endl;
        ber_total.push_back((double)ber/NbChnl);
        //print into a file
        data_out << (SNR_u+9) << ",\t" << (double)ber/NbChnl << endl;
    }
    
    data_out.close();
    /******************************** End of Polar Decoding ***********************/

    //Convert bits to char
    int abits[8];
    //int conv_err=0;

    vector<vector<char> > coded_packets;//(PC_data.size(), vector<char> (PC_data[0].size(), 0));
    for (unsigned int i = 0; i < PC_data.size(); ++i){ //data_bits.size()
        int packetSize = N/8;
        vector<char> conv(packetSize,0);
        //conv_err=0;
        //cout << PC_data[i].size() << ", ";
        for (int j = 0; j < packetSize; ++j){
            for (int k = 0; k < 8; ++k){
                abits[k] = bits_coded[i][j*8+k];//data_bits[i][j*8+k];
                conv[j] += abits[k]*pow(2,k);
            }
        }
        coded_packets.push_back(conv);
        //cout << coded_packets.at(i).size() << " ";
    }
    //cout <<  endl << endl;
    int dataSize = data_bits.size();
    for (int k = 0; k < dataSize; ++k){
        PC_data[k].clear();
    }

    //This function still misses adding the headers
    return coded_packets;

}

PC initialize_PC(int N, int K_w){ //, double SNR_w, double SNR_s

      PC PC_w;
      double designSNRdb = 0;

      // Weak user
      PC_w.constructPC(N, K_w, designSNRdb);

      return PC_w;
}

vector< vector<char> > codingVarCodeRate(int *coloring, int n_col, data_matrix data, cf_data outputForColoring,
    header_transmission **header_data, vector<int> coderate){

    //std::vector<int> chunks_Node;
    string pathFiles = "../repository";//"/CachingFile/repository"

    unsigned int b_chuncks = data.b_chunks;

    nodo *nodi = outputForColoring.nodes;
    int n = outputForColoring.n_nodi;

    //Calculate the number of chunks per node based on the code rate
    //6 because we assume the chunk size is 1/6 of the packet size
    //for (unsigned int i = 0; i < coderate.size(); ++i)
    //    chunks_Node.push_back(round(coderate[i]*6)); NO NEED SINCE THE CODE RATE IS UPDATED TO BE CHUNKS PER NODE

    //open a file -- added new
    //ofstream data_out;
    //data_out.open("Results.txt",ios::app);

    header_transmission *header = new header_transmission[n_col];
    if (!header){
        cout << endl << "MEMORY EXCEPTION: header memory allocation." << endl << endl;
        exit(0);
    }

    vector<string> files;
    files = getDirectoryFiles(pathFiles);

    if (files.empty()){
        cout << endl << "Repository folder is empty or it not exist." << endl << endl;
        exit(0);
    }

    /*Sort file by name*/
    sort(files.begin(), files.end());

    /*Find max size on each package*/
    double max_size = -INF;
    for (string f : files){
        int size_file = getFileSize(f);
        //cout << size_file << ", ";
        if (size_file > max_size)
            max_size = size_file;
    }


    int max_size_package = ceil(max_size / (double) b_chuncks);

    vector< vector<char> > coded_data(n_col, vector<char>(3*max_size_package, 0));

    /*For each transmission*/
    for (int i = 0; i < n_col; i++){
        int color = i + 1;
        int begin_data = 1;
        vector<char> v;

        /*For each node*/
        for (int j = 0; j < n; j++){
            if (coloring[j] == color){

                nodo *single_node = &(nodi)[j];

                unsigned int id_utente = single_node->id_utente;
                unsigned int id_file = single_node->id_file;
                vector<int> id_chunck = single_node->id_chunck;

                cout << "[" << color << ", ";
                for(int kk=0; kk < id_chunck.size() ; kk++)
                    cout << id_chunck[kk] << ", ";
                cout << id_utente << "], ";
                //unsigned int size_pkg;

                /*Reading file */
                string file = files.at(id_file);
                ifstream is (file, ifstream::binary);

                if (is) {
                    int begin_package, size_package;
                    double size_file = getFileSize(file);
                    int size_chunk = ceil(size_file / (double) b_chuncks);
                    //counter introduced on the packet that consists of several concatenated chunks
                    int cc=0;
                    int chunks_size = id_chunck.size();

                    for (int m = 0; m < chunks_size; ++m){
                        begin_package = (id_chunck[m] * size_chunk);

                        if (begin_package >= size_file){
                            begin_package = size_file;
                            size_package = 0;
                        }else{
                            if (begin_package + size_package > size_file){
                                size_package = abs(size_file - begin_package);
                            }else
                                size_package = size_chunk;
                        }
                        //cout << "[" << begin_package << ", " << id_chunck[m] << "], ";

                        if (size_package > 0){
                            is.seekg(begin_package);

                            char *buffer = new char[size_package];

                            is.read(buffer, size_package);



                            for (int k = 0; k < size_package; k++){
                                if (begin_data){
                                    coded_data.at(i).at(cc) = buffer[k];
                                    //if(i<5)
                                    //    cout << coded_data.at(i).at(cc);
                                    //v.at(cc) = buffer[k];
                                }else{
                                    coded_data.at(i).at(cc) = (char) coded_data.at(i).at(cc) ^ buffer[k];
                                    //v.at(cc) = (char) v.at(cc) ^ buffer[k];
                                }
                                cc++;
                            }
                            //if(i<5)
                            //    cout << endl << endl;

                            delete[] buffer;

                        }
                        header[i].id_chunks.push_back(id_chunck[m]);

                        //THE BELOW NOT NECESSARY SINCE THE USER KNOWS THE NUMBER OF CHUNKS IN A FILE
                        // AND THE SIZE OF THE FILE WHICH MEANS CAN EXTRACT THE SIZE OF THE LAST CHUNK
                        //THE FLAG should contain the position of the last chunk value between 1 & 3
                        //If the value of the flag is 0, then there is no last chunk
                        //ADD THE SIZE OF THE CHUNK THAT IS LAST IN CASE FLAG != 0
                        //AT THE RECEPTION, THE SIZE CHUNKS (WHICH CORRESPOND TO LAST CHUNKS)
                        //SHOULD BE READ WISELY, BECAUSE WE ONLY INCLUDE THE SIZE OF THE
                        //CHUNKS THAT ARE LAST, SO WE NEED TO KNOW EACH SIZE CORREPSONDS TO WHICH PACKET
                        //IN THE MIXED ONES

                    }
                    //if (begin_data)
                    //    coded_data.push_back(v);

                    begin_data = 0;
                    header[i].id_utenti.push_back(id_utente);
                    header[i].id_files.push_back(id_file);
                    header[i].nb_chunks.push_back(id_chunck.size());
                    //header[i].lastChunk_flag = flag;
                    //header[i].size_package.push_back(size_package); // NO NEED

                }else{
                        cout << endl << "Error opening file in storage: " << file << endl;
                        exit(0);
                }
                is.close();
            }
        }
    }
    cout << "---Coding for weak user OK---" << endl;

    //print into a file
    /*for (int i = 0; i < n_col; i++)
        data_out << header[i].id_utenti.size() << endl;*/
    //data_out.close();

    (*header_data) = header;
    return coded_data;
}

bool myfunction (int i,int j) {
    return (i<j);
}

vector< vector<char> > codingVarCodeRate_Ref(int *coloring, int n_col, data_matrix data, cf_data outputForColoring,
    header_transmission **header_data, vector<int> coderate){

    //std::vector<int> chunks_Node;
    string pathFiles = "../repository";

    unsigned int b_chuncks = data.b_chunks;

    nodo *nodi = outputForColoring.nodes;
    int n = outputForColoring.n_nodi;
    int u, lim;

    vector<vector <nodo> > nodes_list;//all node lists for each color
    vector<nodo> node_color; //node list for one color
    vector<vector<int> > cr_color;//vector that contains the users code rates in each nodes' list of one color
    vector<vector<int> > cr_cumcolor;//vector that contains the cumulative users code rates in each list of one color
    vector<int> a;
    vector<int> b;
    //two loops to detect the gather the nodes for each color
    for (int i = 1; i <= n_col; ++i)
    {
        a =  vector<int> ();
        b =  vector<int> ();
        node_color = vector<nodo> ();
        int cc = 0, l=0;
        for (int j = 0; j < n; ++j)
        {
            if(i == coloring[j]){
                node_color.push_back(nodi[j]);
                u = nodi[j].id_utente;
                a.push_back(coderate[u]); //add the code rates of the users of each color
                /*for(int ii=0; ii< nodi[j].id_chunck.size(); ii++)
                    cout << nodi[j].id_chunck[ii] << ",";
                cout << endl;*/
                cc++;
            }
        } 
        if(cc > 0){
            sort(a.begin(), a.end(),myfunction);
            b.push_back(0);

            for(unsigned int k=0; k < a.size(); k++){
                if(a[k] != a[k+1]){
                    b.push_back(b[l]+a[k]); //contains the distinct code rates sorted e.g. [0 1 3]
                    l++;
                }
            }
            /*if(i<41){
                cout <<"[";
                for(unsigned int k=0; k < a.size(); k++)
                    cout << a[k] << ",";
                cout <<"], ";
                for(unsigned int k=0; k < b.size(); k++)
                    cout << b[k] << ",";
                cout << endl;
            }*/
            cr_color.push_back(a);
            cr_cumcolor.push_back(b);
            nodes_list.push_back(node_color);
            /*if(i<41)
                for (unsigned int jj = 0; jj < nodes_list[i-1].size(); ++jj)
                    cout << nodes_list[i-1][jj].id_chunck.size() << ", " << nodes_list[i-1].size() << endl;*/
        }

    }

    vector<string> files;
    files = getDirectoryFiles(pathFiles);

    if (files.empty()){
        cout << endl << "Repository folder is empty or it not exist." << endl << endl;
        exit(0);
    }

    /*Sort file by name*/
    sort(files.begin(), files.end());

    /*Find max size on each package*/
    double max_size = -INF;
    for (string f : files){
        int size_file = getFileSize(f);
        //cout << size_file << ", ";
        if (size_file > max_size)
            max_size = size_file;
    }


    int max_size_package = ceil(max_size / (double) b_chuncks);

    /*For each transmission*/
    unsigned int id_utente, id_file, ct=-1;
    int nodes_len = nodes_list.size();
    vector< vector<char> > coded_data;
    vector<header_transmission> header_t;

    for (int kk = 0; kk < nodes_len; ++kk){

        a = cr_color.at(kk);
        b = cr_cumcolor.at(kk);
        //lim is the limit where we should stop when we run over the cumulative code rate vec b
        lim = 0; 
        while(b[lim++] <= a[a.size()-1] && lim < b.size()){ 
        }
        //cout << "[" << lim << ", " << a[a.size()-1] << "], ";

        //loop on the different parts of the code rate that can be go together
        for (int  i = 1; i < lim; ++i){
            //begin_data=1 means that a new subchunk to be introduced, =0 means a subchunk to XOR
            //pre_pak is the size of the previous packet to know where to start the new subchunk
            int begin_data=1, pre_pak=0;
            int begin_package, size_package=0;
            vector<char> v;
            header_transmission h_tmp_u;
            vector<header_transmission> h_tmp (nodes_list[kk].size(), h_tmp_u);//each element is one XORed subchunk
            ct++; //counter on the coded_data
            //cout << ct <<",";

            //loop on the subchunks
            for (int j = b[i-1]; j < b[i]; ++j){
                pre_pak = size_package;
                begin_data=1;
                //loop on the chunks to XOR
                for (unsigned int jj = 0; jj < nodes_list[kk].size(); ++jj){

                    if(nodes_list[kk][jj].id_chunck.size() > j){
                        //subchunk[j] ^=  nodes_list[kk][jj].id_chunks[j];
                        nodo single_node = nodes_list[kk][jj];

                        if(j==b[i-1]){
                            id_utente = single_node.id_utente;
                            id_file = single_node.id_file;
                        }
                        int id_chunck = single_node.id_chunck[j]; //current chunk to be read
                        //cout << "[" << id_chunck << ",";

                        /*Reading file */
                        string file = files.at(id_file);
                        ifstream is (file, ifstream::binary);
                        if (is) {

                            double size_file = getFileSize(file);
                            int size_chunk = ceil(size_file / (double) b_chuncks);

                            begin_package = (id_chunck * size_chunk);

                            if (begin_package >= size_file){
                                begin_package = size_file;
                                size_package = 0;
                             }else{
                                if (begin_package + size_package > size_file){
                                  size_package = abs(size_file - begin_package);
                                }else
                                  size_package = size_chunk;
                            }
                            //cout << "[" << begin_package << ", " << id_chunck << "], ";

                            if (size_package > 0){
                                is.seekg(begin_package);
                                char *buffer = new char[size_package];
                                is.read(buffer, size_package);

                                //cout << ct << ", " << begin_data << "," << j*pre_pak << "," << j << "," << jj << "], ";
                                if(begin_data){
                                    for (int k = 0; k < size_package; k++){
                                        v.push_back(buffer[k]);
                                    }
                                    begin_data = 0;
                                    //cout << v.size() << ", ";
                                }else{
                                    for (int k = 0; k < size_package; k++){
                                        v[j*pre_pak+k] = v[j*pre_pak+k]^buffer[k];
                                    }
                                }
                                delete[] buffer;
                                h_tmp[jj].id_chunks.push_back(id_chunck);
                                if(j==b[i-1]){
                                    h_tmp_u.id_utenti.push_back(id_utente);
                                    h_tmp_u.id_files.push_back(id_file);
                                }
                            }


                        }else{
                            cout << endl << "Error opening file in storage: " << file << endl;
                            exit(0);
                        }
                        is.close();
                    }
                }
            }
            coded_data.push_back(v);
            //cout << coded_data[ct].size() << ", ";
            //cout << nodes_list[kk].size() << ", ";
            for (unsigned int jj = 0; jj < nodes_list[kk].size(); ++jj){
                //cout << h_tmp[jj].id_chunks.size() << ", ";
                if(h_tmp[jj].id_chunks.size() > 0)
                    h_tmp_u.nb_chunks.push_back(h_tmp[jj].id_chunks.size());
                for (int j = 0; j < h_tmp[jj].id_chunks.size(); ++j){
                    h_tmp_u.id_chunks.push_back(h_tmp[jj].id_chunks[j]);
                    //cout << h_tmp[jj].id_chunks[j] << ",";
                }
            }
            /*for(unsigned int jj = 0; jj < h_tmp_u.id_chunks.size(); ++jj)
                cout << h_tmp_u.id_chunks[jj] << ",";*/

            header_t.push_back(h_tmp_u);
            /*for (unsigned int jj = 0; jj < header_t[ct].id_chunks.size(); ++jj)
                cout << header_t[ct].id_chunks[jj] << ", ";
            cout << endl;*/

        }

    }
    //cout << endl  << header_t.size() << ", " << nodes_list.size() << endl;
    /*for(int i=0; i< header_t.size(); i++){
        //int index = index_find(header_t[i].id_utenti, 2);
        //if(index!=-1){
            for(int k=0; k< header_t[i].id_utenti.size(); k++)
                cout << "[" << header_t[i].id_utenti[k] << "," << header_t[i].id_files[k] << "," << header_t[i].nb_chunks[k] << ",";
            for(int k=0; k<header_t[i].id_chunks.size();k++)
                cout << header_t[i].id_chunks[k] << ",";
            cout << "]" << endl;
        //}
    }*/

    cout << "---Coding for weak user OK---" << endl;

    //Copy header vector
    int header_len = header_t.size();
    header_transmission *header = new header_transmission[header_len];
    if (!header){
        cout << endl << "MEMORY EXCEPTION: header memory allocation." << endl << endl;
        exit(0);
    }
    for(int i=0; i < header_len; i++)
        header[i] = header_t.at(i);
    (*header_data) = header;

    return coded_data;
}


void decodingData(header_transmission header, vector<char> &coded_data, unsigned int m_files, unsigned int b_chuncks,
	unsigned int id_user, unsigned int id_demand, unsigned int id_requested_file, int c_rate, unsigned int *n_package_remains){

    string pathFolder = "../cache/UserCache/user_" + to_string(id_user); //"/CachingFile/cache/UserCache/user_"

    char *coded_file_buffer = NULL;
    int trovato = 0;
	int size_chunk = coded_data.size()/c_rate;
	int size_packet = coded_data.size();
	vector<vector<int> > id_chunks;
	unsigned int nb_chunks;
	unsigned int TotalChunks=0;
    int cc;

    //Reconstruct the chunks matrix [users  x chunk ID]
    for (unsigned int i = 0; i < header.id_utenti.size(); i++){
        /*All info for requested package*/
        nb_chunks = header.nb_chunks.at(i);
        vector<int> vd;
		for(unsigned int j=TotalChunks; j< nb_chunks + TotalChunks; j++)
			vd.push_back(header.id_chunks.at(j));
		id_chunks.push_back(vd);
		TotalChunks += nb_chunks;

	}
	//Find the user index
	int index = index_find(header.id_utenti, id_user);
    //cout << index << " ";

	coded_file_buffer = new char[size_packet];
    for (int k = 0; k < size_packet; k++){
        coded_file_buffer[k] = coded_data.at(k);
    }

	for (unsigned int i = 0; i < header.id_utenti.size(); i++){
		if(id_user != header.id_utenti.at(i)){
            //cout << "[" <<  header.id_utenti.at(i) << "," << header.id_files.at(i) << ",";

			/********************************** READ CACHE FILE *********************************/
            char *cache_chunk = new char[size_chunk];
			vector<char> cache (size_packet, 0);
			cc=0;

            int max_cache = (id_chunks.at(i).size()>c_rate)? c_rate:id_chunks.at(i).size();
			for(unsigned int j=0; j < max_cache; j++){ //id_chunks.at(i).size()
                //cout << id_chunks.at(i).at(j) << ",";
				// Open file cache for read package
				string pathFileCache = pathFolder + "/" + to_string(header.id_files.at(i)) + "_" + to_string(id_chunks.at(i).at(j)) + ".cache";

				ifstream inFilePackage (pathFileCache, ifstream::binary);
				if (inFilePackage){
					inFilePackage.read(cache_chunk, size_chunk);
					inFilePackage.close();
				} else{
					cout << endl << "Error reading file cache: " << pathFileCache << endl;
					// exit(0);
				}
				for(int k=0; k<size_chunk; k++){
					cache[cc] = cache_chunk[k];
					cc ++;
				}

			}
            //cout << "]";
            /***********************************************************************************/

			/*Reading packages in cache*/
            for (int k = 0; k < size_packet; k++){
                coded_file_buffer[k] = coded_file_buffer[k] ^ cache[k];
            }
			delete[] cache_chunk;
		}
		if(index != -1)
			trovato = 1;
	}

     /*   if (id_user == id_utente && id_file == id_requested_file){
            (*n_package_remains)--;

            if (size_package_ > 0){

                trovato = 1;
                id_requested_chunck = id_chunck;
                size_requested_package = size_package_;

                coded_file_buffer = new char[coded_data.size()];
                for (unsigned int k = 0; k < coded_data.size(); k++){
                    coded_file_buffer[k] = coded_data.at(k);
                }

                for (unsigned int j = 0; j < header.id_files.size(); j++){
                    unsigned int id_chunck_xor = header.id_chunks.at(j);
                    unsigned int id_file_xor = header.id_files.at(j);
                    unsigned int size_package_xor = header.size_package.at(j);

                    /*************************************************/
					//This case is when the same chunk is requested by several user, so it is sent once for all
					//THE BELOW IS COMMENTED BECAUSE IT IS NOT TAKEN INTO ACCOUNT ON THE SENDER
                    /*int isAlreadyAnalized = 0;
                    unsigned int r = 0;
                    while (r < j && !isAlreadyAnalized){
                        if (id_file_xor == header.id_files.at(r) && id_chunck_xor == header.id_chunks.at(r)){
                            isAlreadyAnalized = 1;
                        }
                        r++;
                    }*/
                    /*************************************************/
                   /* if (i != j && (id_file != id_file_xor || id_chunck != id_chunck_xor) && !isAlreadyAnalized && size_package_xor > 0){
                        /********************************** READ CACHE FILE *********************************/
                    /*    char *cache = new char[size_package_xor];

                        // Open file cache for read package
                        string pathFileCache = pathFolder + "/" + to_string(id_file_xor) + "_" + to_string(id_chunck_xor) + ".cache";

                        ifstream inFilePackage (pathFileCache, ifstream::binary);
                        if (inFilePackage){
                            inFilePackage.read(cache, size_package_xor);
                            inFilePackage.close();
                        }
                        else{
                            cout << endl << "Error reading file cache: " << pathFileCache << endl;
                            exit(0);
                        }
                        /***********************************************************************************/

                        /*Reading packages in cache*/
                    /*    for (unsigned int k = 0; k < size_package_xor; k++){
                            coded_file_buffer[k] = coded_file_buffer[k] ^ cache[k];
                        }
                    }
                }
            }
        }
    }*/

    if (trovato){
        /********************************** DELIVERY PACKAGE *********************************/
        cc=0;
		for(unsigned int i=0; i<id_chunks.at(index).size(); i++){
			// Open file for write delivery package
			string pathFileDelivery = pathFolder + "/" + to_string(header.id_files.at(index)) + "_" + to_string(id_chunks.at(index).at(i)) + ".cache";

			char *chunk = new char[size_chunk];
			for(int j=0; j<size_chunk; j++){
				chunk[j] = coded_file_buffer[cc];
				cc ++;
			}

			ofstream outFileDelivery (pathFileDelivery, ifstream::binary);
			if (outFileDelivery.is_open()){
				outFileDelivery.write(chunk, size_chunk);
				outFileDelivery.close();
				(*n_package_remains)--;
			}
			else{
				cout << endl << "Error writing delivery package: " << pathFileDelivery << endl;
				exit(0);
			}
            delete[] chunk;
		}

        /************************************************************************************/
    }

    //For testing pourpose
    if ((*n_package_remains) < 5 ){
        string name_file;
        name_file = "../trasmissioni/User_" + to_string(id_user) + "/decoded_file_" + to_string(id_demand) + ".xml"; //"/CachingFile/trasmissioni/User_"
        ofstream outFile (name_file, ios::out | ios::binary);

        for (unsigned int i = 0; i < b_chuncks; i++){
            /********************************** READ ALL PACKAGE ********************************/
            // Open file cache for read package
            string pathFilePackage = pathFolder + "/" + to_string(id_requested_file) + "_" + to_string(i) + ".cache";
            int size_file = getFileSize(pathFilePackage);
            ifstream outFilePackage (pathFilePackage, ifstream::binary);

            if (size_file > 0){
                char *buffer = new char[size_file];

                if (outFilePackage){
                    outFilePackage.read(buffer, size_file);
                    outFilePackage.close();
                }
                else{
                    cout << endl << "Error reading package: " << pathFilePackage << endl;
                    // exit(0);
                }
                /**********************************************************************************/

                outFile.write(buffer, size_file);

                delete[] buffer;
            }
        }

        outFile.close();
    }

    delete[] coded_file_buffer;
}

void dyn_chnl(gr_complex out_ch[], gr_complex in_ch[], gr_complex h[], int Nc, int L, int t[]){
    
    //signal with cyclic pefix
    gr_complex in_ch1[t[L-1]+Nc];
    for (int i = 0; i < t[L-1]; ++i){
        in_ch1[i] = in_ch[Nc-t[L-1]+i];

    }

    memcpy((void *) &in_ch1[t[L-1]], (void *) &in_ch[0], sizeof(gr_complex)*Nc);
    for (int i = t[L-1]; i < Nc+t[L-1]; ++i)
    {
        out_ch[i-t[L-1]] = 0;
        for (int l = 0; l < L; ++l)
        {
            out_ch[i-t[L-1]] += h[l]*in_ch1[i-t[l]];
        }
        //out_ch[i-t[L-1]] /= norm;
    }                         
    
}

//======================================================================
gr_complex** generate_ch_param(double pow_del[], int Nb_Symb, int Nc,int Ls,int fd,int T,int L){
    
    double delta_0 = pow(10,-6);
    double delta = delta_0;
    double u;
    double SQRT_Ls = sqrt(Ls);
    double cos_f,sin_f;
    gr_complex sos;
    gr_complex **a;
    double phi[L][Ls];
    double ksi[L][Ls];
    double theta[L][Nb_Symb];
    a = new gr_complex*[Nb_Symb];
    for (int i = 0; i < Nb_Symb; ++i)
        a[i] = new gr_complex[L];

    for (int m = 0; m < L; ++m){
        //cout << "phi and ksi: " ;
        for (int l = 0; l < Ls; ++l)
        {
            phi[m][l] = -pi + 2*((double) rand()/(double)(RAND_MAX))*pi;
            ksi[m][l] = -pi + 2*((double) rand()/(double)(RAND_MAX))*pi;
            //cout << "[" << phi[m][l] << ", " << ksi[m][l] << "], "; 
        }

        theta[m][0] = -pi + ((double) rand()/(double)(RAND_MAX))*2*pi;
        //cout << "theta: " << theta[m][0] << endl;

        for (int j=1; j<=Nb_Symb; ++j){
            u   = ((double) rand()/(double)(RAND_MAX))*pi;           
            if(theta[m][j-1]>=pi) delta = -1*delta_0;
            if(theta[m][j-1]<=pi) delta = delta_0;
            theta[m][j] =  theta[m][j-1] + delta*u;
            theta[m][j] = max((double)-pi, min(theta[m][j],(double)pi));
            //cout << "theta: " << theta[m][j] << endl;
            
            cos_f=0; sin_f=0;
            for (int l = 0; l < Ls; ++l)
            {
                cos_f += cos(2*pi*fd*j*T*cos((double) (2*pi*l-pi+theta[m][j])/(4*Ls) + phi[m][l]));
                sin_f += sin(2*pi*fd*j*T*sin((double) (2*pi*l-pi+theta[m][j])/(4*Ls) + ksi[m][l]));
            }
            //cout << endl;
            sos.real((double) cos_f/SQRT_Ls);
            sos.imag((double) sin_f/SQRT_Ls);
            //cout << "SOS: " << sos << endl;
            a[j-1][m].real(sqrt(pow_del[m])*sos.real());
            a[j-1][m].imag(sqrt(pow_del[m])*sos.imag());
            //cout << a[j-1][m] << ", ";         
            
        }
    }
    return a;
}

//======================================================================

void print( const char * prompt, gr_complex A[], int log2N )
{
   int N = 1 << log2N;
   cout << prompt << '\n' << fixed;
   for ( int i = 0; i < N; i++ ) cout << A[i] << '\n';
}

//======================================================================

void FFT( gr_complex f[], gr_complex ftilde[], int log2N )                 // Fast Fourier Transform
{
   int N = 1 << log2N;

   // Reorder
   for ( int i = 0; i < N; i++ )
   {
      int ii = 0, x = i;
      for (int j = 0; j < log2N; j++)
      {
         ii <<= 1;
         ii |= ( x & 1 );
         x >>= 1;
      }
      ftilde[ii] = f[i];
   }

   for ( int s = 1; s <= log2N; s++ )
   {
      int m = 1 << s;
      int m2 = m / 2;
      gr_complex w = 1.0;
      gr_complex wm = polar( (float) 1.0, -pi / m2 );
      for ( int j = 0; j < m2; j++ )
      {
         for ( int k = j; k < N; k += m )
         {
            gr_complex t = w * ftilde[k+m2];
            gr_complex u =     ftilde[k   ];
            ftilde[k   ] = u + t;
            ftilde[k+m2] = u - t;
         }
         w *= wm;
      }
   }
}

//======================================================================

void iFFT( gr_complex ftilde[], gr_complex f[], int log2N )                // Inverse Fast Fourier Transform
{
   int N = 1 << log2N;

   for ( int m = 0; m < N; m++ ) ftilde[m] = conj( ftilde[m] );      // Apply conjugate (reversed below)

   FFT( ftilde, f, log2N );

   float factor = 1.0 / N;
   for ( int m = 0; m < N; m++ ) f[m] = conj( f[m] ) * factor;

   for ( int m = 0; m < N; m++ ) ftilde[m] = conj( ftilde[m] );      // Only necessary to reinstate ftilde
}
//======================================================================

gr_complex divide_complex(gr_complex s1, gr_complex s2){
    //Division
    double a=(((s1.real())*(s2.real()))+((s1.imag())*(s2.imag())))/(pow(s2.real(),2)+pow(s2.imag(),2));
    double b=(((s2.real())*(s1.imag()))-((s1.real())*(s2.imag())))/(pow(s2.real(),2)+pow(s2.imag(),2));
    gr_complex s;
    s.real(a); s.imag(b);
    return s;

}
//======================================================================


}//end namespace caching
