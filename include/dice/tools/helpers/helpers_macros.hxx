#ifndef __HELPERS_MACROS_HXX__
#define __HELPERS_MACROS_HXX__

// Remove warnings for unused variable
#define UNUSED_VARIABLE(expr) (void)(expr)

#define OFFSETOF(type, field)    ((size_t) ( (char *)&((type *)(0))->field - (char *)0 ))
#define OFFSETOF_D(inst,field)    ((size_t) ( (char *)&inst->field - (char *)&inst ))

#define DO_STRINGIFY(x) #x
#define STRINGIFY(x) DO_STRINGIFY(x)
#define CONCAT_ONE_ARGUMENT(x,arg) x ( arg )
/*
#define MAKE_ENUM(name, ...) enum class name { __VA_ARGS__, __COUNT}; \
inline std::ostream& operator<<(std::ostream& os, name value) { \
std::string enumName = #name; \
std::string str = #__VA_ARGS__; \
int len = str.length(); \
std::vector<std::string> strings; \
std::ostringstream temp; \
for(int i = 0; i < len; i ++) { \
if(isspace(str[i])) continue; \
        else if(str[i] == ',') { \
        strings.push_back(temp.str()); \
        temp.str(std::string());\
        } \
        else temp<< str[i]; \
} \
strings.push_back(temp.str()); \
os << enumName << "::" << strings[static_cast<int>(value)]; \
return os;} 
*/
//#define COMMA() ,

//#define OFFSETOF(type, field)    ((size_t) &(((type *) 0)->field))
//#define OFFSETOF(type, field) offsetof(type,field)

#endif
