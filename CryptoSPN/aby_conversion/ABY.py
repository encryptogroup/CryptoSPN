"""
Created on January 9, 2019

@author: Amos Treiber
"""
import subprocess
import tempfile

from spn.structure.Base import Product, Sum, Leaf, get_nodes_by_type
from spn.structure.leaves.parametric.Parametric import Bernoulli, Poisson, Gaussian
from spn.structure.leaves.histogram.Histograms import Histogram
from math import log, e, pi

from .CPPConversion import aby_header, aby_footer, dec_str, Selection, cmake_file
from CryptoSPN.Constants import CRYPTOSPN_DIR, ABY_DIR, COMPILE_TIMEOUT

import logging

logger = logging.getLogger(__name__)


def spn_to_aby(node, bitlen=64, leaves=None, sel=Selection.OSELECTION, aby_inputs="", aby_circuit="",
               aby_declarations=None, nodes_done=None):
    if nodes_done is None:
        nodes_done = set()
    if aby_declarations is None:
        aby_declarations = {}
    if leaves is None:
        leaves = {}
    if node in nodes_done:
        return bitlen, leaves, aby_inputs, aby_circuit, aby_declarations, nodes_done

    if isinstance(node, Bernoulli):
        pstr = "std::numeric_limits<" + dec_str(bitlen) + ">::lowest()" if node.p == 0 else str(log(node.p, 2))
        opstr = "std::numeric_limits<" + dec_str(bitlen) + ">::lowest()" if node.p == 1 else str(log(1 - node.p, 2))

        aby_inputs += "        p = " + pstr + ";\n        op = " + opstr + ";\n"
        aby_inputs += "        pptr = (uint" + str(bitlen) + "_t*) &p;\n" + "        opptr = (uint" + str(
            bitlen) + "_t*) &op;\n"
        aby_inputs += "        share* s_p_" + str(node) + " = circ->PutINGate(*pptr, bitlen, SERVER);\n" \
                                                          "        share* s_op_" + str(node) \
                      + " = circ->PutINGate(*opptr, bitlen, SERVER);\n"

        aby_declarations.update({"p": (dec_str(bitlen), ";\n"), "op": (dec_str(bitlen), ";\n"),
                                 "pptr": ("uint" + str(bitlen) + "_t*", ";\n"),
                                 "opptr": ("uint" + str(bitlen) + "_t*", ";\n")})

        if sel == Selection.LINEAR:
            aby_inputs += "        idx = " + str(node.scope[0]) + ";\n"
            aby_inputs += "        share* s_idx_" + str(
                node) + " = circ->PutINGate(idx, 32, SERVER,);\n"
            aby_circuit += "        share* s_" + str(node) + " = putBernoulliNode(putInputSelector(s_idx_" + str(
                node) + ", s_V_in, s_ids, circ), s_op_" + str(node) + ", s_p_" + str(node) + ", circ);\n"
            aby_declarations.update({"s_V_in": ("std::vector<share*>", """;
        for (uint64_t b : v_V) {
            s_V_in.push_back(circ->PutINGate(b, 1, CLIENT));
        }\n"""), "s_ids": ("std::vector<share*>", """;
        for (auto i=0; i < v_V.size(); i++) {
            s_ids.push_back(circ->PutCONSGate((uint" + str(bitlen) + "_t) i, 32));
        }

        uint32_t idx;\n""")})
        elif sel == Selection.NONE:
            aby_circuit += "        share* s_" + str(node) + " = putBernoulliNode(s_V_in.at(" + str(
                node.scope[0]) + "), s_op_" + str(node) + ", s_p_" + str(node) + ", circ);\n"
            aby_declarations.update({"s_V_in": ("std::vector<share*>", """;
        for (uint64_t b : v_V) {
            s_V_in.push_back(circ->PutINGate(b, 1, CLIENT));
        }\n""")})
        elif sel == Selection.OSELECTION:
            aby_circuit += "        share* s_" + str(node) + " = putBernoulliNode(sel[" + str(
                leaves[node]) + "], s_op_" + str(node) + ", s_p_" + str(node) + ", circ);\n"
        else:
            raise Exception("Selection Type not recognized")
        nodes_done.add(node)
        return bitlen, leaves, aby_inputs, aby_circuit, aby_declarations, nodes_done

    if isinstance(node, Gaussian):
        log2e = log(e, 2)

        aby_inputs += "        mu = " + str(node.mean) + ";\n"
        aby_inputs += "        log2ps2 = " + str(0.0 - log(2.0 * pi * node.variance, 2)) + ";\n"
        aby_inputs += "        loge2s2 = " + str(log2e / (2.0 * node.variance)) + ";\n"
        aby_inputs += "        muptr = (uint" + str(bitlen) + "_t*) &mu;\n" + "        log2ps2ptr = (uint" + \
                      str(bitlen) + "_t*) &log2ps2;\n" + "        loge2s2ptr = (uint" + str(bitlen) + "_t*) &loge2s2;\n"
        aby_inputs += "        share* s_mu_" + str(node) + " = circ->PutINGate(*muptr, bitlen, SERVER);\n" \
                                                           "        share* s_log2ps2_" + \
                      str(node) + " = circ->PutINGate(*log2ps2ptr, bitlen, SERVER);\n     " \
                                  "   share* s_loge2s2_" + str(node) + \
                      " = circ->PutINGate(*loge2s2ptr, bitlen, SERVER);\n"

        aby_declarations.update(
            {"mu": (dec_str(bitlen), ";\n"), "log2ps2": (dec_str(bitlen), ";\n"), "loge2s2": (dec_str(bitlen), ";\n"),
             "muptr": ("uint" + str(bitlen) + "_t*", ";\n"), "log2ps2ptr": ("uint" + str(bitlen) + "_t*", ";\n"),
             "loge2s2ptr": ("uint" + str(bitlen) + "_t*", ";\n")})

        if sel == Selection.LINEAR:
            aby_inputs += "        idx = " + str(node.scope[0]) + ";\n"
            aby_inputs += "        share* s_idx_" + str(node) + \
                          " = circ->PutINGate(idx, 32, SERVER);\n"
            aby_circuit += "        share* s_" + str(node) + " = putLogGaussianNode(putInputSelector(s_idx_" + \
                           str(node) + ", s_V_in, s_ids, circ), s_mu_" + str(node) + ", s_log2ps2_" + str(node) + \
                           ", s_loge2s2_" + str(node) + ", circ);\n"
            aby_declarations.update({"s_ids": ("std::vector<share*>", """;
        for (auto i=0; i < v_V.size(); i++) {
            s_ids.push_back(circ->PutCONSGate((uint64_t) i, 32));
        }

        uint32_t idx;\n"""), "s_V_in": ("std::vector<share*>", """;
        """ + dec_str(bitlen) + """ val;
        for (double k : v_V) {
            val = k;
            s_V_in.push_back(circ->PutINGate(*(uint""" + str(bitlen) + """_t*) &val, bitlen, CLIENT));
        }\n""")})

        elif sel == Selection.NONE:
            aby_circuit += "        share* s_" + str(node) + " = putLogGaussianNode(s_V_in.at(" + str(
                node.scope[0]) + "), s_mu_" + str(node) + ", s_log2ps2_" + str(node) + ", s_loge2s2_" + str(
                node) + ", circ);\n"
            aby_declarations.update({"s_V_in": ("std::vector<share*>", """;
        """ + dec_str(bitlen) + """ val;
        for (double k : v_V) {
            val = k;
            s_V_in.push_back(circ->PutINGate(*(uint""" + str(bitlen) + """_t*) &val, bitlen, CLIENT));
        }\n""")})

        elif sel == Selection.OSELECTION:
            aby_circuit += "        share* s_" + str(node) + " = putLogGaussianNode(sel[" + str(
                leaves[node]) + "], s_mu_" + str(node) + ", s_log2ps2_" + str(node) + ", s_loge2s2_" + str(
                node) + ", circ);\n"

        else:
            raise Exception("Selection Type not recognized")

        nodes_done.add(node)
        return bitlen, leaves, aby_inputs, aby_circuit, aby_declarations, nodes_done

    if isinstance(node, Poisson):
        log2e = log(e, 2)
        llstr = "std::numeric_limits<" + dec_str(bitlen) + ">::min()" if node.mean == 0 else str(log(node.mean, 2))

        aby_inputs += "        loglambda = " + llstr + ";\n"
        aby_inputs += "        lambdaloge = " + str(node.mean * log2e) + ";\n"
        aby_inputs += "        llptr = (uint" + str(bitlen) + "_t*) &loglambda;\n" + "        lleptr = (uint" + \
                      str(bitlen) + "_t*) &lambdaloge;\n"
        aby_inputs += "        share* s_ll_" + str(node) + \
                      " = circ->PutINGate(*llptr, bitlen, SERVER);\n        share* s_lle_" + \
                      str(node) + " = circ->PutINGate(*lleptr, bitlen, SERVER);\n"

        aby_declarations.update({"loglambda": (dec_str(bitlen), ";\n"), "lambdaloge": (dec_str(bitlen), ";\n"),
                                 "llptr": ("uint" + str(bitlen) + "_t*", ";\n"),
                                 "lleptr": ("uint" + str(bitlen) + "_t*", ";\n")})

        if sel == Selection.LINEAR:
            aby_inputs += "        idx = " + str(node.scope[0]) + ";\n"
            aby_inputs += "        share* s_idx_" + str(node) + \
                          " = circ->PutINGate(idx, 32, SERVER, seed_exp, &expander);\n"
            aby_circuit += "        share* s_" + str(node) + " = putLogPoissonNode(putInputSelector(s_idx_" + \
                           str(node) + ", s_V_in, s_ids, circ), s_ll_" + str(node) + ", putInputSelector(s_idx_" + \
                           str(node) + ", s_V_logfac, s_ids, circ), s_lle_" + str(node) + ", circ);\n"
            aby_declarations.update({"s_ids": ("std::vector<share*>", """;
        for (auto i=0; i < v_V.size(); i++) {
            s_ids.push_back(circ->PutCONSGate((uint64_t) i, 32));
        }

        uint32_t idx;\n"""), "s_V_in": ("std::vector<share*>", """;
        std::vector<share*> s_V_logfac;
        """ + dec_str(bitlen) + """ val, val_logfac;
        for (double k : v_V) {
            val = k;
            val_logfac = log2(tgamma(k + 1));
            s_V_logfac.push_back(circ->PutINGate(*(uint""" + str(bitlen) + """_t*) &val_logfac, bitlen, CLIENT));
            s_V_in.push_back(circ->PutINGate(*(uint""" + str(bitlen) + """_t*) &val, bitlen CLIENT));
        }\n""")})

        elif sel == Selection.NONE:
            aby_circuit += "        share* s_" + str(node) + " = putLogPoissonNode(s_V_in.at(" + str(
                node.scope[0]) + "), s_ll_" + str(node) + ", s_V_logfac.at(" + str(node.scope[0]) + "), s_lle_" + str(
                node) + ", circ);\n"
            aby_declarations.update({"s_V_in": ("std::vector<share*>", """;
        std::vector<share*> s_V_logfac;
        """ + dec_str(bitlen) + """ val, val_logfac;
        for (double k : v_V) {
            val = k;
            val_logfac = log2(tgamma(k + 1));
            s_V_logfac.push_back(circ->PutINGate(*(uint""" + str(bitlen) + """_t*) &val_logfac, bitlen, CLIENT));
            s_V_in.push_back(circ->PutINGate(*(uint""" + str(bitlen) + """_t*) &val, bitlen, CLIENT));
        }\n""")})

        elif sel == Selection.OSELECTION:
            aby_circuit += "        share* s_" + str(node) + " = putLogPoissonNode(sel[" + str(
                leaves[node]) + "], s_ll_" + str(node) + ", sel_logfac[" + str(node.scope[0]) + "], s_lle_" + str(
                node) + ", circ);\n"
            aby_declarations.update({"v_V_logfac": ("std::vector<double>", """;
        """ + dec_str(bitlen) + """ val;
        for (double k : v_V) {
            v_V_logfac.push_back(log2(tgamma(k + 1)));
        }\n""")})

        else:
            raise Exception("Selection Type not recognized")

        nodes_done.add(node)
        return bitlen, leaves, aby_inputs, aby_circuit, aby_declarations, nodes_done

    if isinstance(node, Histogram):
        aby_declarations.update(
            {"v_s_borders": ("std::vector<share*>", ";\n		for (double k : {" + ", ".join(map(str, node.breaks))
                             + """}) {
            """ + dec_str(bitlen) + """ val = k;
            v_s_borders.push_back(circ->PutINGate(*(uint""" + str(bitlen) +
                             """_t*) &val, bitlen, SERVER));
        }\n""")})
        aby_inputs += "        std::vector<" + dec_str(bitlen) + "> v_" + str(node) + "_densities = {" + ", ".join(
            map(str, map(lambda d: log(d, 2), node.densities))) + "};\n"
        aby_inputs += "        std::vector<uint" + str(bitlen) + "_t> v_" + str(node) + "_densitiesi = dtoi(v_" + str(
            node) + "_densities);\n"
        aby_inputs += "        std::vector<share*> v_s_" + str(node) + "_densities;\n"
        aby_inputs += "        for(uint" + str(bitlen) + "_t k : v_" + str(node) + """_densitiesi) {
            v_s_""" + str(node) \
                      + """_densities.push_back(circ->PutINGate(k, bitlen, SERVER));
        }\n"""

        if sel == Selection.NONE:
            aby_declarations.update({"s_V_in": ("std::vector<share*>", """;
        for (double in : v_V) {
            """ + dec_str(bitlen) + """ val = in;
            s_V_in.push_back(circ->PutINGate(*(uint""" + str(bitlen)
                                                + """_t*) &val, bitlen, CLIENT));
        }\n""")})
            aby_circuit += "        share* s_" + str(node) + " = putHistogramNode(s_V_in.at(" + \
                           str(node.scope[0]) + "), v_s_borders, v_s_" + str(node) + "_densities, circ);\n"

        elif sel == Selection.LINEAR:
            aby_inputs += "        idx = " + str(node.scope[0]) + ";\n"
            aby_inputs += "        share* s_idx_" + str(
                node) + " = circ->PutINGate(idx, 32, SERVER);\n"
            aby_circuit += "        share* s_" + str(node) + " = putHistogramNode(putInputSelector(s_idx_" + str(
                node) + ", s_V_in, s_ids, circ), v_s_borders, v_s_" + str(node) + "_densities, circ);\n"
            aby_declarations.update({"s_V_in": ("std::vector<share*>", """;
        for (uint64_t in : v_V) {
            """ + dec_str(bitlen) + """ val = in;
            s_V_in.push_back(circ->PutINGate(val, bitlen, CLIENT));
        }\n"""), "s_ids": ("std::vector<share*>", """;
        for (auto i=0; i < v_V.size(); i++) {
            s_ids.push_back(circ->PutCONSGate((uint32_t) i, 32));
        }

        uint32_t idx;\n""")})

        elif sel == Selection.OSELECTION:
            aby_circuit += "        share* s_" + str(node) + " = putHistogramNode(sel[" + str(leaves[node]) + \
                           "], v_s_borders, v_s_" + str(node) + "_densities, circ);\n"

        else:
            raise Exception("Selection Type not recognized")

        nodes_done.add(node)
        return bitlen, leaves, aby_inputs, aby_circuit, aby_declarations, nodes_done

    if isinstance(node, Product) or isinstance(node, Sum):
        res_list = list(map(lambda child: spn_to_aby(child, bitlen, leaves, sel, aby_inputs, aby_circuit,
                                                     aby_declarations, nodes_done), node.children))
        aby_inputs = "".join(map(lambda aby: aby[2], res_list))
        aby_circuit = "".join(map(lambda aby: aby[3], res_list))
        map(lambda aby: aby_declarations.update(aby[4]), res_list)
        map(lambda aby: nodes_done.add(aby[5]), res_list)

        if isinstance(node, Product):
            aby_circuit += "        std::vector<share*> v_" + str(node) + "_children = { " + ", ".join(
                map(lambda child: "s_" + str(child), node.children)) + " };\n"
            aby_circuit += "        share* s_" + str(node) + " = putLogProdNode(v_" + str(node) + "_children, circ);\n"

        else:
            aby_inputs += "        std::vector<" + dec_str(bitlen) + "> v_" + str(node) + "_weights = { " + ", ".join(
                map(lambda w: str(log(w, 2)), node.weights)) + " };\n"
            aby_inputs += "        std::vector<uint" + str(bitlen) + "_t> v_" + str(
                node) + "_weightsi = dtoi( v_" + str(
                node) + "_weights);\n"
            aby_inputs += "        std::vector<share*> v_" + str(node) + "_sweights;\n"
            aby_inputs += "        for (int j = 0; j < v_" + str(node) + """_weights.size(); ++j) {
            v_""" + str(node) + "_sweights.push_back(circ->PutINGate(v_" + str(node) \
                          + """_weightsi.at(j), bitlen, SERVER));
        }\n"""

            aby_circuit += "        std::vector<share*> v_" + str(node) + "_children = { " + ", ".join(
                map(lambda child: "s_" + str(child), node.children)) + " };\n"
            aby_circuit += "        share* s_" + str(node) + " = putLogSumNode(v_" + str(node) + "_children, v_" + str(
                node) + "_sweights, v_" + str(node) + "_children.size() , circ);\n"

        nodes_done.add(node)
        return bitlen, leaves, aby_inputs, aby_circuit, aby_declarations, nodes_done

    raise Exception("Node type not registered: " + str(type(node)))


def spn_to_aby_file(node, cryptospn_path=CRYPTOSPN_DIR, bitlen=64, filename="spntest.cpp", sel=Selection.OSELECTION):
    logger.info(f"Creating {filename}...")
    aby_head = aby_header(cryptospn_path, bitlen)

    aby_end = aby_footer(node, bitlen, filename)

    if sel == Selection.OSELECTION:
        (bitlen, leaves, aby_inputs, aby_circuit, aby_declarations, nodes_done) = spn_to_aby(node, bitlen,
                                                                                             {node: pos for pos, node in
                                                                                              enumerate(
                                                                                                  get_nodes_by_type(
                                                                                                      node, Leaf))},
                                                                                             sel)
        selection_input = "        std::vector<uint32_t>ids = {"
        for node, pos in leaves.items():
            selection_input += str(node.scope[0]) + ", "

        if isinstance(next(iter(leaves.keys())), Poisson):
            selection_input += "};\n        std::vector<uint" + str(
                bitlen) + "_t > v_ConvV = dtoi(v_V);\n        std::vector<uint" + str(
                bitlen) + "_t > v_ConvV_logfac = dtoi(v_V_logfac);\n" \
                          "        share** sel= selection_GC(v_ConvV, ids.size(), ids.data(), circ);\n" \
                          "        share** sel_logfac = selection_GC(v_ConvV_logfac, ids.size(), ids.data(), circ);\n"

        elif isinstance(next(iter(leaves.keys())), Gaussian) or isinstance(next(iter(leaves.keys())), Histogram):
            selection_input += "};\n        std::vector<uint" + str(bitlen) \
                               + "_t > v_ConvV = dtoi(v_V);\n" \
                                 "        share** sel= selection_GC(v_ConvV, ids.size(), ids.data(), circ);\n "

        else:
            selection_input += "};\n        std::vector<uint" + str(bitlen) \
                               + "_t> v_ConvV(v_V.begin(), v_V.end());\n" \
                                 "        share** sel= selection_GC(v_ConvV, ids.size(), ids.data(), circ);\n"
        aby_inputs = selection_input + aby_inputs

    else:
        (bitlen, _, aby_inputs, aby_circuit, aby_declarations, nodes_done) = spn_to_aby(node, bitlen, sel=sel)

    for key, (before, after) in aby_declarations.items():
        aby_inputs = "        " + before + " " + key + after + aby_inputs

    f = open(filename, 'w')

    f.write(aby_head + aby_inputs + aby_circuit + aby_end)
    logger.info(f"{filename} created.")


def spn_to_aby_exec(node, aby_path=ABY_DIR, cryptospn_path=CRYPTOSPN_DIR, bitlen=64, name="spntest",
                    sel=Selection.OSELECTION):
    spn_to_aby_file(node, cryptospn_path, bitlen, name + ".cpp", sel)

    with tempfile.NamedTemporaryFile() as tmpfile:
        cmake_filename = tmpfile.name
        with open(cmake_filename, 'w') as f:
            f.write(cmake_file(name + ".cpp", aby_path, name))

        logger.debug(f"CMakeLists.txt created in {cmake_filename}.")

        logger.info(f"Compiling {name}. This might take some time...")
        if cryptospn_path[-1] != '/':
            cryptospn_path += '/'
        proc = subprocess.Popen([cryptospn_path + 'compiling/compile.sh', cryptospn_path, aby_path, name + '.cpp',
                                 name, cmake_filename])

        ret = proc.wait(timeout=COMPILE_TIMEOUT)
        if ret == 0:
            logger.info(f"Created executable {name}.")
        elif ret == 2:
            logger.warning(f"Compilation of {name} failed! Did you forget to build ABY?")
        else:
            logger.warning(f"Compilation of {name} failed!")
