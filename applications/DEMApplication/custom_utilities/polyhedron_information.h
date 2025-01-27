/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

#if !defined(POLYHEDRON_INFORMATION_H)
#define  POLYHEDRON_INFORMATION_H

namespace Kratos
{

class PolyhedronInformation
{
    
public:
    PolyhedronInformation(){};
    std::string mName;
    std::vector<double> mListOfSize;
    std::vector<double> mListOfVolume;
    std::vector<std::vector<array_1d<double,3>>> mListOfVerticesList;
    std::vector<std::vector<std::vector<int>>> mListOfFacesList;
    std::vector<array_1d<double,3>> mListOfInertiaPerUnitMass;
    
    virtual ~PolyhedronInformation() {};
    
    virtual void PrintInfo(std::ostream& rOStream) const {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {
        //rOStream << std::endl << this->mSize <<"  " << this->mVolume << std::endl;                                
    }
    
    virtual std::string Info() const {
        std::stringstream buffer;
        buffer << mName;
        return buffer.str();
    }  
    
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
         //TODO:   
    }

    void load(Serializer& rSerializer)
    {   
        //TODO:            
    }		

};

inline std::istream& operator >> (std::istream& rIStream, PolyhedronInformation& rThis);

inline std::ostream& operator << (std::ostream& rOStream, const PolyhedronInformation& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#endif // POLYHEDRON_INFORMATION_H  defined 