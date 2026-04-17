/* -*- c++ -*- */


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "CheckFunction.h"
#include "EnvironmentSetup.h"
#include "DecodingInput.h"
#include "ConflictGraph.h"
#include "randomHandler.h"
#include "DataDefinition.h"
#include "grasp.h"
#include "CodingDecodingData.h"
#include "FuncsFromMain.h"
#include "TxRx.h"
//#include "TestTxRx.h"
#include <iomanip>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <random>

using namespace caching;
using namespace std;

//const float pi = 3.1415926;

/*void print( const char * prompt, gr_complex A[], int log2N );
void FFT ( gr_complex f[]     , gr_complex ftilde[], int log2N );
void iFFT( gr_complex ftilde[], gr_complex f[]     , int log2N );
gr_complex** generate_ch_param(double pow_del[], int OFDM_Symb, int Nc,int Ls,int fd,int T,int L);
void dyn_chnl(gr_complex out_ch[], gr_complex in_ch[], gr_complex h[], int Nc, int L, int t[]);
*/
int main()
{
    int m_files = 20;
    int b_chunks = 200;
    int id_demand = 0;
    int nb_strg = 1;
    int id_user = 0;
    int size_chunk =42;
    data_matrix d_data;
    cf_data d_outputForColoring;
    bool d_GRASP=true;
    int d_n_col=0;
    int *d_coloring;
    //bool DEBUG = true;
    PC PC_w, PC_s;
    float designSNRdb = 0;
    unsigned int packet_remain;
    bool isStr;
    int packet_len = 256;

    header_transmission *d_header_data;
    header_transmission *d_hdr_sdata;

    vector<int> d_spack_size;
    vector< vector<char> > d_coded_data;
    vector< vector<char> > d_strg_data;
    vector< vector<char> > d_PC_data;
    vector<header_polar> d_hX; //this is the main header to be added directly to the transmitted packets
    vector<vector<char> > d_transmission;
    vector<int> header_len;
    //vector<unsigned int> d_small_packet_size;
    std::vector<int> coderate = {1, 2, 2, 3, 4};
    int K_w = coderate[id_user]*size_chunk*8;
    int K_s = 4*size_chunk*8;
    int N=2048;

    //srand((unsigned) time(NULL));

    cout << endl << "Data generation process" << endl << "-------------" << endl << endl;

    d_data = generateData(m_files, b_chunks, id_demand);
    int nb_users = d_data.n_utenti;
    /*for(int i=0; i<nb_users; i++){
        for(int j=0; j<d_b_chunks; j++)
            cout << d_data.Q_chuncks[i][j]  << ", ";
        cout << endl << endl;
    } */


    cout << endl << "Conflict-Graph generator process" << endl << "-------------" << endl << endl;
    d_outputForColoring = conflictGraphGenerator(d_data, coderate);
    cout << endl << "Numero nodi del grafo: " << d_outputForColoring.n_nodi << endl << endl;

    if (d_outputForColoring.n_nodi > 0)
    {
        cout << endl << "Graph Coloring process" << endl << "-------------" << endl << endl;
        // Coloring
        if (d_GRASP)
        {
            int d_maxIter = 20;
            d_coloring = graspGraphColoring(d_maxIter, d_outputForColoring, &d_n_col);
            colorRienumeration(d_n_col, &d_coloring, d_outputForColoring.n_nodi);
        }

        cout << endl << "La colorazione e' stata effettuata con successo!" << endl;
        cout << endl << "Numero di colori utilizzati: " << d_n_col << endl;
        cout << "The expected gain is: " << (100*(d_outputForColoring.n_nodi-d_n_col)/d_outputForColoring.n_nodi) << "%" << endl;

    }/*end if (d_outputForColoring.n_nodi > 0)*/

    if (d_n_col > 0)
    {
        cout << endl << "Coding data process" << endl << "-------------" << endl << endl;
        //Coding data to be transmitted for weak users - xor coding
        d_coded_data = codingVarCodeRate_Ref(d_coloring, d_n_col, d_data, d_outputForColoring, &d_header_data,coderate);
        //d_coded_data = codingVarCodeRate(d_coloring, d_n_col, d_data, d_outputForColoring, &d_header_data,coderate);
        //d_coded_data = codingData(d_coloring, d_n_col, d_data, d_outputForColoring, &d_header_data);

        //Build the schema (graph) packets for strong and weak users
        vector<vector<bool> > G_edges;
        //bool **G_edges;
        d_strg_data = MaxBipartiteGraph(d_coloring, d_n_col, d_outputForColoring.nodes, 
            d_outputForColoring.n_nodi, nb_strg, d_data, &d_hdr_sdata, G_edges);
        
        G_edges = vector<vector<bool> > ();

        for(unsigned int i=0; i<d_coded_data.size(); i++){
            vector<bool> temp(d_strg_data.size(), false);
            G_edges.push_back(temp);
        }

        //Coding polarly the weak and strong packets
        vector<vector<int> >  bits_coded;
        d_PC_data = codingDataPolar(d_coded_data, d_strg_data, bits_coded, G_edges, d_header_data, d_hdr_sdata, d_hX, N);
        cout << "The total number of transmitted packet is: " << d_PC_data.size() << endl;

        //Pack the data and attach the header to be ready for transmission
        /*TX_simul(d_hX, d_PC_data, d_transmission, header_len);
        
        int TXdata_size = d_transmission.size();
        //vector<int> data_bits (TXdata_size*8, 0);
        vector<vector<int> > data_bits;

        //=============Convert a stream of chars to stream of bits=================/
        std::vector<unsigned int> bb (8,0);

        for(int i=0; i<TXdata_size; i++){
            int packetSize = d_transmission[i].size();
            data_bits.push_back(vector<int>(8*packetSize, 0));

            for(int j=0; j<packetSize; j++){
                bb = conv_char_to_bitsInt(d_transmission[i][j]);
                for(int k=0; k<8; k++){
                    data_bits[i][j*8+k] = bb[k];
                }
            }
            int sum=0;
            if(i<20){
                for (int j = header_len[i]*8; j < data_bits[i].size(); ++j)
                    sum += data_bits[i][j];
                cout << sum << ", ";
            }
        }
        cout << endl;

        //============QPSK mapping=============================
        vector<vector<gr_complex> > tx_Symb;
        tx_Symb = BitsToQPSKSymb(data_bits);
        vector<vector<gr_complex> > rx_Symb;

        for(unsigned int i=0; i<header_len.size();i++)
            header_len[i] *= 4; 

        for (int i = 0; i < tx_Symb.size(); ++i){
            vector<gr_complex> v;
            for (int j = 0; j < tx_Symb[i].size(); ++j)
                v.push_back(tx_Symb[i][j]);
            rx_Symb.push_back(v);
        }


        


        //============Polar codes environment=============================
        PC_w = initialize_PC(N,K_w);
        PC_s.initPC(N, K_s, designSNRdb);
        PC_s.setGenMatrix(PC_w.genMatrix);
        PC_s.setRn(PC_w.Rn);
        PC_s.setArrangedBits(PC_w.arragedBits);

        //==============================================================
        //=---------Apply the IFFT block for OFDM multiplexing---------=
        //=---------Channel block                             ---------=
        //=---------Apply the FFT block for OFDM demultiplexing--------=
        //==============================================================
        const int N_carr = 64, log2N = 6;       // Hardwired for testing
        gr_complex h       [N_carr];//Channel coefficients in Frequency domain 
        gr_complex f       [N_carr];//original symbols in Frequency domain 
        gr_complex yf      [N_carr];//received symbols in Frequency domain 
        gr_complex in_ch   [N_carr];//original symbols in Time domain
        gr_complex out_ch  [N_carr];//received symbols in Time domain
        gr_complex out_ch_p  [N_carr];//noiseless received symbols in Time domain
        int OFDM_Symb = ceil(packet_len*4/N_carr);
        //gr_complex a[L][N_carr]; for all users

        int Nb_packets = tx_Symb.size();
        //generate channel elements independently of j time slot
        int L=3, Ls=8, fd=40, T=80*pow(10,-6); 
        int t[] = {0, 2, 4};
        gr_complex **a;//[Nb_packets][L];//fo 1 user, to be modified
        double pow_delay[] = {0.5, 0.3, 0.2};
        gr_complex noise;
        
        int chnl_test=1, snrmin = 20, snrmax=20;
        double sqrtVariance, snr;
        vector<double> ber ((snrmax-snrmin+1), 0);
        int cctt=0;

        for(snr = snrmin; snr <= snrmax; snr++){
            sqrtVariance = sqrt(pow(10,-(snr/10)));
            for(int ch=0; ch<chnl_test; ch++)
            {

                a = generate_ch_param(pow_delay, OFDM_Symb*Nb_packets, N_carr, Ls, fd, T, L);

                //noise generator function
                default_random_engine generator( rand() );
                normal_distribution<double> distribution(0,sqrtVariance);

                
                for (int i = 0; i < Nb_packets; ++i)
                {
                    for (int j = 0; j < OFDM_Symb; ++j)
                    {
                        // Forward transform
                        memcpy((void *) &f[0], (void *) &tx_Symb[i][header_len[i]+j*N_carr], sizeof(gr_complex)*N_carr);

                        //Invert FFT for OFDM multiplexing
                        iFFT( f, in_ch, log2N );          

                        //Channel block
                        //generate attenuation vector a1,a2,...,aL for one OFDM symb i.e. length N_carr
                        dyn_chnl(out_ch_p, in_ch, a[i*OFDM_Symb+j], N_carr, L, t);

                        //PC_w.noise_gc(out_ch_p, out_ch, sqrtVariance,N_carr);//

                        //FFT of noiseless received signal for channel estimation
                        if(j==0){
                            FFT ( out_ch_p, yf, log2N );   //print( "\nTransform:", f, log2N );
                            //Channel estimation and equalization
                            //cout << "Channel coefficients: ";
                            float tot_norm =0;
                            for (int k = 0; k < N_carr; ++k){
                                h[k]= yf[k]/f[k];
                                tot_norm += pow(abs(h[k]),2);
                            }
                            tot_norm /= N_carr;
                            for (int k = 0; k < N_carr; ++k)
                                h[k] /= sqrt(tot_norm);

                        }
                        
                        //FFT for OFDM demultiplexing
                        FFT ( out_ch_p, yf, log2N );   //print( "\nTransform:", f, log2N );
                        //Equalization
                        //cout << "\nEqualized signal: ";
                        for (int k = 0; k < N_carr; ++k){
                            yf[k]= yf[k]/h[k];
                        }


                        for(int k=0; k<N_carr; ++k){
                            noise = gr_complex(distribution(generator)/sqrt(2),distribution(generator)/sqrt(2));
                            yf[k] = yf[k];// + noise/h[k];
                        }

                        //SREBUILD THE RECEIVED DATA SYMBOLS FROM yt
                        memcpy((void *) &rx_Symb[i][header_len[i]+j*N_carr], (void *) &yf[0], sizeof(gr_complex)*N_carr);
                        
                    }
                    memcpy((void *) &rx_Symb[i][0], (void *) &tx_Symb[i][0], sizeof(gr_complex)*header_len[i]);
                    
                }
                
                
                //============Decoding environment=============================
                
                if(id_user < d_data.n_utenti)
                    packet_remain = readCacheInfo(id_user, 0);
                else 
                    packet_remain = b_chunks;
                
                if ((d_data.n_utenti-1) < id_user)
                    isStr = true;
                else
                    isStr = false;
                
                int ct_err = 0;
                double Pe_In_Pck;
                double Pe=0;
                int nbPacket = 5;//rx_Symb.size();
                for(int i=0; i<nbPacket; i++){
                    Pe_In_Pck = -1;
                    RX_simul(rx_Symb[i], packet_remain, coderate, id_user, PC_w, PC_s, b_chunks, m_files, size_chunk, N, snr, isStr,i,bits_coded[i],Pe_In_Pck);
                    if(Pe_In_Pck != -1){
                        Pe += Pe_In_Pck;
                        ct_err++;
                    }
                }
                Pe /= ct_err;
                cout << Pe << ", ";
                ber[cctt] += Pe;

            }
            ber[cctt] /= chnl_test;
            cout << "Average Error probability: " << ber[cctt] << endl;
            cctt++;

        }
        cout << "Final BER: " << endl;
        for(int i=0; i< ber.size(); i++)
            cout << ber[cctt] << ", ";
        cout << endl << endl;

        int len_a = OFDM_Symb*Nb_packets;
        for (int i = 0; i < len_a; ++i)
            delete [] a[i];

    */}/* end if (d_n_col > 0) */
    
    delete [] d_coloring;
        
    return 0;

}
/*

void dyn_chnl(gr_complex out_ch[], gr_complex in_ch[], gr_complex h[], int Nc, int L, int t[]){
    
    //signal with cyclic pefix
    gr_complex in_ch1[t[L-1]+Nc];
    for (int i = 0; i < t[L-1]; ++i){
        in_ch1[i] = in_ch[Nc-t[L-1]+i];

    }
    /*cout << "Channel time domain: ";
    for (int l = 0; l < L; ++l)
        cout << pow(abs(h[l]),2) << ", ";
    cout << endl;*/
    /*float norm = 0;
    for (int l = 0; l < L; ++l)
        norm += pow(abs(h[l]),2);
    norm = sqrt(norm);
    //cout << norm << ", ";*/

/*    memcpy((void *) &in_ch1[t[L-1]], (void *) &in_ch[0], sizeof(gr_complex)*Nc);
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
*/