#ifndef DEBUG_H
#define DEBUG_H
    #define DEBUG2
    #if defined DEBUG4 
    #   define debug4(x) std::cout<<"[TRACE] "<<x<<std::endl;
    #   define debug3(x) std::cout<<"[DEBUG] "<<x<<std::endl;
    #   define debug2(x) std::cout<<"[INFO] "<<x<<std::endl;
    #   define debug1(x) std::cout<<"[WARNING] "<<x<<std::endl;
    #   define debug0(x) std::cout<<"[ERROR] "<<x<<std::endl;
    #elif defined DEBUG3 
    #   define debug4(x) 
    #   define debug3(x) std::cout<<"[DEBUG] "<<x<<std::endl;
    #   define debug2(x) std::cout<<"[INFO] "<<x<<std::endl;
    #   define debug1(x) std::cout<<"[WARNING] "<<x<<std::endl;
    #   define debug0(x) std::cout<<"[ERROR] "<<x<<std::endl;
    #elif  defined DEBUG2 
    #   define debug4(x) 
    #   define debug3(x) 
    #   define debug2(x) std::cout<<"[INFO] "<<x<<std::endl;
    #   define debug1(x) std::cout<<"[WARNING] "<<x<<std::endl;
    #   define debug0(x) std::cout<<"[ERROR] "<<x<<std::endl;
    #elif defined DEBUG1 
    #   define debug4(x) 
    #   define debug3(x) 
    #   define debug2(x)
    #   define debug1(x) std::cout<<"[WARNING] "<<x<<std::endl;
    #   define debug0(x) std::cout<<"[ERROR] "<<x<<std::endl;
    #elif defined DEBUG0
    #   define debug4(x) 
    #   define debug3(x) 
    #   define debug2(x)
    #   define debug1(x) 
    #   define debug0(x) std::cout<<"[ERROR] "<<x<<std::endl;
    #else
    #   define debug4(x)
    #   define debug3(x)
    #   define debug2(x)
    #   define debug1(x)
    #   define debug0(x)
    #endif
#endif
