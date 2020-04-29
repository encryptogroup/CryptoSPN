"""
@author: Amos Treiber
"""
from enum import Enum


class Selection(Enum):
    NONE = 0
    LINEAR = 1
    OSELECTION = 2


def dec_str(bitlen):
    if bitlen == 32:
        return "float"
    else:
        return "double"


def aby_header(cryptospn_dir, bitlen):
    return """
#include <ENCRYPTO_utils/crypto/crypto.h>
#include <ENCRYPTO_utils/parse_options.h>
#include "../../abycore/aby/abyparty.h"
#include "../../abycore/circuit/share.h"
#include "../../abycore/circuit/booleancircuits.h"
#include "../../abycore/sharing/sharing.h"
#include "selection_blocks/e_SelectionBlock.h"
#include "selection_blocks/t_SelectionBlock.h"
#include <cassert>
#include <iomanip>
#include <iostream>
#include <random>
#include <limits>
#include <math.h>

void read_test_options(int32_t* argcp, char*** argvp, e_role* role,
                       uint32_t* bitlen, uint32_t* nvals, uint32_t* secparam, std::string* address,
                       uint16_t* port, int32_t* test_op, uint32_t* sharing, uint32_t* iter, std::string* input_file) {

    uint32_t int_role = 0, int_port = 0, int_testbit = 0;

    parsing_ctx options[] =
            { {(void*) &int_role, T_NUM, "r", "Role: 0/1", true, false },
              {(void*) nvals, T_NUM, "n",	"Number of parallel operation elements", false, false },
              {(void*) bitlen, T_NUM, "b", "Bit-length, default 32", false,false },
              {(void*) secparam, T_NUM, "s", "Symmetric Security Bits, default: 128", false, false },
              {(void*) address, T_STR, "a", "IP-address, default: localhost", false, false },
              {(void*) &int_port, T_NUM, "p", "Port, default: 7766", false, false },
              {(void*) test_op, T_NUM, "t", "Single test (leave out for all operations), default: off", false, false },
              {(void*) sharing, T_NUM, "y", "Type of Sharing used, 0: S_BOOL, 1: S_YAO, default: 1", false, false },
              {(void*) iter, T_NUM, "i",	"Number of iterations for evaluation, default: 1", false, false },
              {(void*) input_file, T_STR, "f",	"Input file of user RVs, default: all_data.txt", false, false },
            };

    if (!parse_options(argcp, argvp, options,
                       sizeof(options) / sizeof(parsing_ctx))) {
        print_usage(*argvp[0], options, sizeof(options) / sizeof(parsing_ctx));
        std::cout << "Exiting" << std::endl;
        exit(0);
    }

    assert(int_role < 2);
    *role = (e_role) int_role;

    if (int_port != 0) {
        assert(int_port < 1 << (sizeof(uint16_t) * 8));
        *port = (uint16_t) int_port;
    }

}

std::vector<uint""" + str(bitlen) + """_t> dtoi(std::vector<float>& in) {
    std::vector<uint""" + str(bitlen) + """_t> res;
    for(auto const& d : in) {
        uint""" + str(bitlen) + "_t* di = (uint" + str(bitlen) + """_t*) &d;
        res.push_back(*di);
    }
    return res;
}

std::vector<uint""" + str(bitlen) + """_t> dtoi(std::vector<double>& in) {
    std::vector<uint""" + str(bitlen) + """_t> res;
    for(auto const& d : in) {
        """ + dec_str(bitlen) + """ df = d;
        uint""" + str(bitlen) + "_t* di = (uint" + str(bitlen) + """_t*) &df;
        res.push_back(*di);
    }
    return res;
}

share* toSIMD(std::vector<share*>& s, uint32_t num, BooleanCircuit* bc) {
    uint32_t i,j;

    uint32_t len = """ + str(bitlen) + """;
    std::vector<std::vector<uint32_t> > res(len);
    std::vector<uint32_t> v;

    for (i = 0; i < num; ++i)  {
        for (j = 0; j < len; ++j) {
            res[j].push_back(s[i]->get_wires()[j]);
        }
    }

    for (j = 0; j < len; ++j) {
        v.push_back(bc->PutCombinerGate(res[j]));
    }

    return new boolshare(v, bc);
}

share* putFPSumGate(share* s_simd, uint32_t num, BooleanCircuit* bc) {
    std::vector<share*>s(num);
    share * s_res;
    for (uint32_t i = 0; i < num; ++i) {
        uint32_t pos[] = {i};
        s_res = i == 0 ? ((Circuit*) bc)->PutSubsetGate(s_simd, pos, 1) :
                bc->PutFPGate(s_res, ((Circuit*) bc)->PutSubsetGate(s_simd, pos, 1), ADD, """ + str(bitlen) \
           + """, 1, no_status);
    }
    return s_res;
}

share* putProdNode(std::vector<share*>& in, BooleanCircuit* circ) {
    double one = 1.0;
    uint64_t *oneptr = (uint64_t*) &one;
    share* s_res = circ->PutCONSGate(*oneptr, """ + str(bitlen) + """);

    for(auto const& s_in: in) {
        s_res = circ->PutFPGate(s_res, s_in, MUL, """ + str(bitlen) + """, 1, no_status);
    }

    return s_res;
}

share* putLogProdNode(std::vector<share*>& v_in, BooleanCircuit* circ) {
    return putFPSumGate(toSIMD(v_in, v_in.size(), circ), v_in.size(), circ);
}

share* putSumNode(std::vector<share*> v_in, share* s_simdw, uint32_t num, BooleanCircuit* circ) {
    share* s_res= circ->PutFPGate(toSIMD(v_in, num, circ), s_simdw, MUL, """ + str(bitlen) + """, num, no_status);

    s_res = putFPSumGate(s_res, num, circ);

    return s_res;
}

share* putLogSumNodeSIMD(std::vector<share*> v_in, share* s_simdw, uint32_t num, BooleanCircuit* circ) {
    share* s_res= circ->PutFPGate(toSIMD(v_in, num, circ), s_simdw, ADD, """ + str(bitlen) + """, num, no_status);

    s_res = circ->PutFPGate(s_res, EXP2, """ + str(bitlen) + """, num, no_status);

    s_res = putFPSumGate(s_res, num, circ);

    s_res = circ->PutFPGate(s_res, LOG2, """ + str(bitlen) + """, 1, no_status);

    return s_res;
}

share* putLogSumNode(std::vector<share*> v_in, std::vector<share*> v_w, uint32_t num, BooleanCircuit* circ) {
    share* s_res;
    if (num == 1) {
        s_res = circ->PutFPGate(v_in.at(0), v_w.at(0), ADD, """ + str(bitlen) + """, 1, no_status);
    } else {
        std::vector<share*> v_prod;

        for (int i = 0; i < v_in.size(); ++i) {
            v_prod.push_back(circ->PutFPGate(circ->PutFPGate(v_in.at(i), v_w.at(i), ADD, """ + str(
        bitlen) + """, 1, no_status), EXP2, """ + str(bitlen) + """, 1, no_status));
        }

        s_res = putFPSumGate(toSIMD(v_prod, num, circ), num, circ);
        s_res = circ->PutFPGate(s_res, LOG2, """ + str(bitlen) + """, 1, no_status);
    }

    return s_res;
}

share* putBernoulliNode(share* s_x, share* s_0p, share* s_1p, BooleanCircuit* circ) {
    return circ->PutMUXGate(s_1p, s_0p, s_x);
}

share* putLogPoissonNode(share* s_k, share* s_loglambda, share* s_logkfac, share* s_lambdaloge, BooleanCircuit* circ) {
    share* s_res = circ->PutFPGate(s_k, s_loglambda, MUL, """ + str(bitlen) + """, 1, no_status);

    s_res = circ->PutFPGate(s_res, s_logkfac, SUB, """ + str(bitlen) + """, 1, no_status);
    s_res = circ->PutFPGate(s_res, s_lambdaloge, SUB, """ + str(bitlen) + """, 1, no_status);

    return s_res;
}

share* putLogGaussianNode(share* s_x, share* s_mu, share* s_log2ps2, share* s_loge2s2, BooleanCircuit* circ) {
    share* s_res = circ->PutFPGate(s_x, s_mu, SUB, """ + str(bitlen) + """, 1, no_status);

    s_res = circ->PutFPGate(s_res, SQR, """ + str(bitlen) + """, 1, no_status);

    s_res = circ->PutFPGate(s_res, s_loge2s2, MUL, """ + str(bitlen) + """, 1, no_status);
    s_res = circ->PutFPGate(s_log2ps2, s_res, SUB, """ + str(bitlen) + """, 1, no_status);

    return s_res;
}

share* putHistogramNode(share* s_in, std::vector<share*>& v_borders, std::vector<share*>& v_values, BooleanCircuit* circ) {
    share* s_res = circ->PutMUXGate(v_values.at(0),  v_values.at(1), circ->PutGTGate(v_borders.at(1), s_in));
    for (int i = 2; i < v_values.size(); ++i) {
        s_res = circ->PutMUXGate(s_res, v_values.at(i), circ->PutGTGate(v_borders.at(i), s_in));
    }

    return s_res;
}


share* putInputSelector(share* s_sel, std::vector<share*>& v_in, std::vector<share*>& v_ids, BooleanCircuit* circ) {
    uint""" + str(bitlen) + """_t zero = 0;
    share * s_zero = circ->PutCONSGate(zero, """ + str(bitlen) + """);
    share* s_tmp = circ->PutMUXGate(v_in.at(0), s_zero, circ->PutEQGate(s_sel, v_ids.at(0)));

    for (auto i = 1; i < v_in.size(); i++) {
        s_tmp = circ->PutXORGate(s_tmp, circ->PutMUXGate(v_in.at(i), s_zero, circ->PutEQGate(s_sel, v_ids.at(i))));
    }

    return s_tmp;
}

/**
 * Selection function (garbled circuit)
 * Author: Masoud Naderpour
 */
share** selection_GC(std::vector<uint""" + str(bitlen) + "_t> &featureVec, uint" + str(bitlen) \
           + """_t numDecisionNodes, uint32_t* permutation, BooleanCircuit* &Circ) {

    uint16_t dim = featureVec.size();
    uint16_t m_numNodes = numDecisionNodes;
    uint64_t maxbitlen = """ + str(bitlen) + """;

    //---- init selectionBlcok ---------------
    SelectionBlock *selBlock;
    if (m_numNodes >= dim){
        selBlock = new e_SelectionBlock(dim, m_numNodes, Circ); // extended SelectionBlock
    } else {
        selBlock = new t_SelectionBlock(dim, m_numNodes, Circ); // truncated SelectionBlock
    }

    selBlock->SelectionBlockProgram(permutation);
    selBlock->SetControlBits();

    share **featureVecShr;

    //----------------Setting client input-------------
    featureVecShr = (share**) malloc(sizeof(share*) * dim);
    for(int j = 0; j < dim; j++) {
        featureVecShr[j] = Circ->PutSIMDINGate(1, featureVec[j], maxbitlen, CLIENT);
    }

    //---------------Building selection circuit--------------
    share **out, *cmp, *tmp;

    assert(Circ->GetCircuitType() == C_BOOLEAN);

    std::vector<std::vector<uint32_t> > Inputs(dim);
    for(int i = 0; i < dim; i++){
        Inputs[i].resize(1);
        Inputs[i][0] = Circ->PutCombinerGate(featureVecShr[i]->get_wires());
    }

    std::vector<std::vector<uint32_t> > tempvec(m_numNodes);
    tempvec = selBlock->buildSelectionBlockCircuit(Inputs);

    share **SelectionBlockOutput = (share**) malloc(sizeof(share*) * m_numNodes);
    for (int i = 0; i < tempvec.size(); i++){
        SelectionBlockOutput[i] = new boolshare(Circ->PutSplitterGate(tempvec[i][0]), Circ);
    }

    return SelectionBlockOutput;
}

std::vector<std::vector<double >> readInputs(std::string infilename) {
    std::ifstream infile;
    infile.open(infilename);

    std::string line;
    std::vector<std::vector<double>> res;
    if (infile.is_open()) {
        while (std::getline(infile, line)) {
            if (line.rfind("#", 0) == 0) {
              // skip commented lines
              continue;
            }
            std::vector<double> linevals;
            std::size_t found = 0;
            std::size_t prev_found = found;
            while (found!=std::string::npos) {
                found = line.find(";", found + 1);
                if (found!=std::string::npos) {
                    linevals.push_back(std::stod(line.substr(prev_found, found - prev_found)));
                    prev_found = found + 1;
                }
            }
            //add last element
            std::string lastelem = line.substr(prev_found, line.length() - prev_found);
            if (!lastelem.empty())
                linevals.push_back((std::stod(lastelem)));
            res.push_back(linevals);
        }
    }
    infile.close();

    return std::move(res);
}

//show thousands separator
std::string printBytes(long int  in) {
    std::string s = std::to_string(in);
    for (long int  i = s.length()-3; i>0 ;i=i-3)	{
        s.insert(i, ",");
    }
    return s;
}

void test_spn(e_role role, const std::string& address, uint16_t port, seclvl seclvl, uint32_t nvals, uint32_t nthreads, 
    e_mt_gen_alg mt_alg, e_sharing sharing, uint32_t iter, std::vector<std::vector<double>>& values,
    std::string outfilename) {

    std::ofstream outfile;
    outfile.open(outfilename);

    uint32_t bitlen = """ + str(bitlen) + """;

    ABYParty* party = new ABYParty(role, address, port, seclvl, bitlen, nthreads, mt_alg, 65536 , \"""" + \
           cryptospn_dir + """aby_files/circ");

    double online_time=0, total_time=0;
    uint64_t online_bytes_sent=0, total_bytes_sent=0, online_bytes_recv=0, total_bytes_recv=0;

    for (unsigned int i = 0 ; i < values.size() && i < iter; i++) {
        std::vector<double> v_V = values.at(i);

        std::vector<Sharing*>& sharings = party->GetSharings();


        BooleanCircuit* circ = (BooleanCircuit*) sharings[sharing]->GetCircuitBuildRoutine();\n"""


def aby_footer(node, bitlen, filename):
    return """
        // output gate
        share* res_out = circ->PutOUTGate(s_""" + str(node) + """, ALL);

        // run SMPC
        party->ExecCircuit();

        // retrieve plain text output
        uint32_t out_bitlen, out_nvals;
        uint64_t *out_vals;
        res_out->get_clear_value_vec(&out_vals, &out_bitlen, &out_nvals);

        // print every output
        for (uint32_t i = 0; i < nvals; i++) {

            // dereference output value as double without casting the content
            """ + dec_str(bitlen) + " value = *((" + dec_str(bitlen) + """*) &out_vals[i]);

            value = value * log(2);

            if (outfile.is_open())
                outfile << std::scientific << std::setprecision (std::numeric_limits<""" + dec_str(bitlen) \
           + """>::digits10 + 3) << value << std::endl;
        }

        online_time += party->GetTiming(P_ONLINE);
        total_time += party->GetTiming(P_TOTAL);
        online_bytes_sent += party->GetSentData(P_ONLINE);
        total_bytes_sent += party->GetSentData(P_TOTAL);
        online_bytes_recv += party->GetReceivedData(P_ONLINE);
        total_bytes_recv += party->GetReceivedData(P_TOTAL);

        party->Reset();
    }

    outfile.close();

    std::string role_str;
    switch(role) {
        case CLIENT:
            role_str = "Client.";
            break;
        case SERVER:
            role_str = "Server.";
            break;
    }

    std::cout << "Evaluated " << iter << " spns as " << role_str << std::endl;
    std::cout << "Mean online time was " << online_time/iter << " ms and mean total time was " << total_time/iter
        << " ms." << std::endl;
    std::cout << "Mean online sent was " << printBytes(online_bytes_sent/iter) << " Bytes and mean total sent was " 
        << printBytes(total_bytes_sent/iter) << " Bytes." << std::endl;
    std::cout << "Mean online received was " << printBytes(online_bytes_recv/iter) << " Bytes and mean total recv was "
        << printBytes(total_bytes_recv/iter) << " Bytes." << std::endl;
    std::cout << "=> Mean online Communication was " << printBytes(online_bytes_sent/iter + online_bytes_recv/iter) 
        << " Bytes, " << online_time/iter << " ms and mean total communication was " << 
        printBytes(total_bytes_sent/iter + total_bytes_recv/iter) << " Bytes, " << total_time/iter << " ms." << 
        std::endl;

}

int main(int argc, char** argv) {
    e_role role;
    uint32_t bitlen = """ + str(bitlen) + """, nvals = 1, secparam = 128, nthreads = 1, sharing=1, iter=1;

    uint16_t port = 7766;
    std::string address = "127.0.0.1";
    int32_t test_op = -1;
    e_mt_gen_alg mt_alg = MT_OT;
    std::string input_file = "all_data.txt";

    read_test_options(&argc, &argv, &role, &bitlen, &nvals, &secparam, &address,
                      &port, &test_op, &sharing, &iter, &input_file);

    seclvl seclvl = get_sec_lvl(secparam);

    std::vector<std::vector<double> > values = readInputs(input_file);

    test_spn(role, address, port, seclvl, nvals, nthreads, mt_alg, (e_sharing) sharing, iter, values, \"""" \
           + filename + """_output_data" + std::to_string(role) + ".txt");

    return 0;
}\n"""


def cmake_file(filename, aby_path, exec_name="spntest"):
    return """
cmake_minimum_required(VERSION 3.12)
project(""" + exec_name + """)

set(CMAKE_CXX_STANDARD 14)

find_package(ABY QUIET)
if(ABY_FOUND)
    message(STATUS "Found ABY")
elseif (NOT ABY_FOUND AND NOT TARGET ABY::aby)
    message("ABY was not found: add ABY subdirectory")
    add_subdirectory(""" + aby_path + """  build)
endif()

add_executable(""" + exec_name + " " + filename \
           + """ selection_blocks/e_SelectionBlock.cpp selection_blocks/permutation_network.cpp """ \
           + """selection_blocks/t_SelectionBlock.cpp)
target_link_libraries(""" + exec_name + """ ABY::aby)

    """
