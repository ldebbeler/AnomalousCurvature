//file: testBubble.h
#ifndef TESTBUBBLE_H
#define TESTBUBBLE_H

#include "bubble.h"
#include "datastructs.h"

class testBubble{
    bubble m_bubble;
public:
    testBubble(bubble foo);
    void writeBubble(scalingValues& scale, bubbleValues& bv);
        

};


#endif //TESTBUBBLE_H
