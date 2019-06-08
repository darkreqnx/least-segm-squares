#ifndef SLS_H
#define SLS_H

#include "basic.h"
#include "graph.h"
#include "main.h"

class SLS
{
public:
    SLS();

    //!  Graham Scan Function
    /*!
      A function that takes no parameters. Completes the Graham Scan and stores iterations into a vector of graphs.
      Returns the count of iterations

      \return integer of the iteration count.
    */
    static int SLS_algo(double cost);


    //!  Graham Scan Iteration
    /*!
      A function that takes one integer parameter and sets the global graph object to the same as the index.
      returns integer to show status of success.
      \param x integer that is the index for the current iteration
      \return integer indicating success
    */
    static int Iteration(int x);
};

#endif // SLS_H
