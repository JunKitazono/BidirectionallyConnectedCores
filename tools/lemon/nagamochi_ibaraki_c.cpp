/*
* 
* this function caluculates mincut of the undirected graph
* the input is adjecent matrix and it returns the charactaristic vector of
* mincut and the mincut value
*  
* this functions uses lemon library
*
* 
* Yuma Aoki, 2020
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

#include <sstream>
#include <math.h>

#include <lemon/smart_graph.h>
#include <lemon/adaptors.h>
// #include <lemon/concepts/digraph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
#include <lemon/nagamochi_ibaraki.h>
#include <lemon/full_graph.h>


using namespace matlab::data;
using matlab::mex::ArgumentList;
using matlab::data::ArrayFactory;
using namespace lemon;
using namespace std;


class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
//         check validity of arguments
        checkArguments(outputs, inputs);
        typedef double Value;
        typedef FullGraph::Node Node;
        ArrayFactory factory;
        
        
        matlab::data::TypedArray<double> input_array = std::move(inputs[0]);
//         get array size
        int size = (int)sqrt((double)input_array.getNumberOfElements());
        
//         create a complete graph
        FullGraph graph(size);
        FullGraph::EdgeMap<double> cap(graph);
        FullGraph::NodeMap<bool> cut(graph);

//         set weight for all arcs
        for (int i = 0; i < size; i ++){
            const Node a = graph.operator()(i);
            for (int j = i+1; j < size; j++){
                if (i != j){
                    const Node b = graph.operator()(j);
                    cap[graph.edge(a,b)] = input_array[i][j];
                }
            }
        }

//         run algorithm
        NagamochiIbaraki<FullGraph, FullGraph::EdgeMap<double>> ni(graph, cap);
        ni.run();
        double minval = (double)ni.minCutMap(cut);
        
//          create charactaristic array of mincut
        matlab::data::TypedArray<int> z_mip = factory.createArray<int>({1, (unsigned long)size});
        for (int i = 0; i < size; i++){
            const Node a = graph.operator()(i);
            if (cut[a]){
                z_mip[i] = 1;
            } else {
                z_mip[i] = 0;
            }
        }
        
//         set outputs
        outputs[0] = z_mip;
        outputs[1] = factory.createScalar<double>(minval);
       
    }
    
//     checking num if square number
    bool checkSquareNumber(int num){
        double db_num = num;
        double sqrt_db_num = sqrt(db_num);
        int sqrt_int_num = sqrt_db_num;
        int int_num = sqrt_int_num * sqrt_int_num;
        return (int_num == num);
    }

// check validity of arguments
    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        if (inputs.size() != 1) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Too many arguments") }));
        }
        if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[0].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(u"error", 
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input matrix must be type double") }));
        }
         if (inputs[0].getDimensions().size() != 2) {
             matlabPtr->feval(u"error", 
                 0, std::vector<matlab::data::Array>({ factory.createScalar("Input must be m-by-n dimension") }));
         }
        if (!checkSquareNumber((int)inputs[0].getNumberOfElements())){
            matlabPtr->feval(u"error", 
                 0, std::vector<matlab::data::Array>({ factory.createScalar("Input must be m-by-m dimension") }));        
        }
    }
};