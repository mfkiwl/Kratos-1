//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#pragma once

// System includes
#include <map>
#include <set>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

// External includes

// Project includes
#include "includes/define.h"
#include "input_output/logger.h"
#include "containers/flags.h"
#include "containers/array_1d.h"
#include "containers/weak_pointer_vector.h"

#define KRATOS_SERIALIZATION_DIRECT_LOAD(type)                          \
    void load(std::string const & rTag, type& rValue)                   \
    {                                                                   \
        load_trace_point(rTag);                                         \
            read(rValue);                                               \
    }                                                                   \
    void load(std::string const & rTag, type const& rValue)             \
    {                                                                   \
        load_trace_point(rTag);                                         \
            read(const_cast<type&>(rValue));                            \
    }                                                                   \
                                                                        \
    void load_base(std::string const & rTag, type& rValue)              \
    {                                                                   \
          load_trace_point(rTag);                                       \
      read(rValue);                                                     \
    }

#define KRATOS_SERIALIZATION_DIRECT_SAVE(type)                          \
    void save(std::string const & rTag, type const & rValue)            \
    {                                                                   \
      save_trace_point(rTag);                                           \
      write(rValue);                                                    \
    }                                                                   \
                                                                        \
    void save_base(std::string const & rTag, type const & rValue)       \
    {                                                                   \
      save_trace_point(rTag);                                           \
      write(rValue);                                                    \
    }

#define KRATOS_SERIALIZATION_DIRECT_CREATE(type)                        \
    void* create(std::string const & rTag, type* prototype)             \
    {                                                                   \
      type* p_new = new type;                                           \
      load(rTag, *p_new);                                               \
      return p_new;                                                     \
    }

#define KRATOS_SERIALIZER_MODE_BINARY                                   \
    if(!mTrace) {
#define KRATOS_SERIALIZER_MODE_ASCII                                    \
    } else {
#define KRATOS_SERIALIZER_MODE_END                                      \
    }
namespace Kratos
{

class ModelPart;
class VariableData;
template <class TDataType> class Variable;


///@name Kratos Classes
///@{

/**
 * @class Serializer
 *
 * \ingroup KratosCore
 *
 * @brief The serialization consists in storing the state of an object into a storage format like data file or memory buffer and also retrieving the object from such a media.
 *
 * @details The serialization consists in storing the state of an object into a storage format like data file or memory buffer and also retrieving the object from such a media.
 * The idea of serialization is based on saving all object's data consecutively in the file or buffer and then load it in the same order.
 * In Kratos a serialization mechanism is used for creating the restart file. So for storing an object into restart file and retrieve it afterward on must add the necessary component used by serialization.
 *
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) Serializer
{
public:
    ///@name  Enum's
    ///@{

    enum PointerType {SP_INVALID_POINTER, SP_BASE_CLASS_POINTER, SP_DERIVED_CLASS_POINTER};
    enum TraceType {SERIALIZER_NO_TRACE=0, SERIALIZER_TRACE_ERROR=1, SERIALIZER_TRACE_ALL=2};

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Serializer
    KRATOS_CLASS_POINTER_DEFINITION(Serializer);

    KRATOS_DEFINE_LOCAL_FLAG( MPI );
    KRATOS_DEFINE_LOCAL_FLAG( SHALLOW_GLOBAL_POINTERS_SERIALIZATION );

    typedef std::size_t SizeType;

    typedef void* (*ObjectFactoryType)();

    typedef std::map<void*, void*> LoadedPointersContainerType;

    typedef std::map<std::string, ObjectFactoryType> RegisteredObjectsContainerType;

    typedef std::map<std::string, std::string> RegisteredObjectsNameContainerType;

    typedef std::set<const void*> SavedPointersContainerType;

    typedef std::iostream BufferType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit Serializer(BufferType* pBuffer, TraceType const& rTrace=SERIALIZER_NO_TRACE) :
        mpBuffer(pBuffer), mTrace(rTrace), mNumberOfLines(0)
    {
    }

    /// Destructor.
    virtual ~Serializer()
    {
        delete mpBuffer;
    }

    ///@}
    ///@name Operators
    ///@{

    /// Sets the Serializer in a state ready to be loaded
    /// Note: If the same object is loaded twice before deleting it from memory all its pointers will be duplicated.
    void SetLoadState();

    ///@}
    ///@name Operations
    ///@{
    ///This function returns the "trace type" used in initializing the serializer.
    ///Trace type is one of SERIALIZER_NO_TRACE,SERIALIZER_TRACE_ERROR,SERIALIZER_TRACE_ALL
    TraceType GetTraceType() const {return mTrace;}

    void SetBuffer(BufferType* pBuffer)
    {
        mpBuffer = pBuffer;
    }

    template<class TDataType>
    static void* Create()
    {
        return new TDataType;
    }

    template<class TDataType>
    static void Register(std::string const & rName, TDataType const& pPrototype)
    {
        msRegisteredObjects.insert(RegisteredObjectsContainerType::value_type(rName,Create<TDataType>));
        msRegisteredObjectsName.insert(RegisteredObjectsNameContainerType::value_type(typeid(TDataType).name(), rName));
    }

    static void Deregister(const std::string& rName)
    {
        msRegisteredObjects.erase(rName);
        msRegisteredObjectsName.erase(rName);
    }

    template<class TDataType>
    void load(std::string const & rTag, TDataType& rObject)
    {
        load_trace_point(rTag);
        rObject.load(*this);
    }

    template<class TDataType>
    void load(std::string const & rTag, Kratos::shared_ptr<TDataType>& pValue)
    {
        PointerType pointer_type = SP_INVALID_POINTER;
        void* p_pointer;
        read(pointer_type);

        if(pointer_type != SP_INVALID_POINTER)
        {
            read(p_pointer);
            LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
            if(i_pointer == mLoadedPointers.end())
            {
                if(pointer_type == SP_BASE_CLASS_POINTER)
                {
                    if constexpr (!std::is_abstract_v<TDataType>) {
                        if(!pValue) {
                            pValue = this->MakeShared<TDataType>();
                        }
                    }
                    else {
                        KRATOS_ERROR << "Cannot instantiate an abstract class\n";
                    }
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    KRATOS_ERROR_IF(i_prototype == msRegisteredObjects.end())
                        << "There is no object registered in Kratos with name : "
                        << object_name << std::endl;

                    if(!pValue) {
                        pValue = Kratos::shared_ptr<TDataType>(static_cast<TDataType*>((i_prototype->second)()));
                    }
                }

                // Load the pointer address before loading the content
                mLoadedPointers[p_pointer]=&pValue;
                load(rTag, *pValue);
            }
            else
            {
                pValue = *static_cast<Kratos::shared_ptr<TDataType>*>((i_pointer->second));
            }
        }
    }

    template<class TDataType>
    void load(std::string const & rTag, Kratos::intrusive_ptr<TDataType>& pValue)
    {
        PointerType pointer_type = SP_INVALID_POINTER;
        void* p_pointer;
        read(pointer_type);

        if(pointer_type != SP_INVALID_POINTER)
        {
            read(p_pointer);
            LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
            if(i_pointer == mLoadedPointers.end())
            {
                if(pointer_type == SP_BASE_CLASS_POINTER)
                {
                    if constexpr (!std::is_abstract_v<TDataType>) {
                        if(!pValue) {
                            pValue = Kratos::intrusive_ptr<TDataType>(new TDataType);
                        }
                    }
                    else {
                        KRATOS_ERROR << "Cannot instantiate an abstract class\n";
                    }
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    KRATOS_ERROR_IF(i_prototype == msRegisteredObjects.end())
                        << "There is no object registered in Kratos with name : "
                        << object_name << std::endl;

                    if(!pValue) {
                        pValue = Kratos::intrusive_ptr<TDataType>(static_cast<TDataType*>((i_prototype->second)()));
                    }
                }

                // Load the pointer address before loading the content
                mLoadedPointers[p_pointer]=&pValue;
                load(rTag, *pValue);
            }
            else
            {
                pValue = *static_cast<Kratos::intrusive_ptr<TDataType>*>((i_pointer->second));
            }
        }
    }

    template<class TDataType>
    void load(std::string const & rTag, Kratos::unique_ptr<TDataType>& pValue)
    {
        PointerType pointer_type = SP_INVALID_POINTER;
        void* p_pointer;
        read(pointer_type);

        if(pointer_type != SP_INVALID_POINTER)
        {
            read(p_pointer);
            LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
            if(i_pointer == mLoadedPointers.end())
            {
                if(pointer_type == SP_BASE_CLASS_POINTER)
                {
                    if constexpr (!std::is_abstract_v<TDataType>) {
                        if(!pValue) {
                            pValue = Kratos::unique_ptr<TDataType>(new TDataType);
                        }
                    }
                    else {
                        KRATOS_ERROR << "Cannot instantiate an abstract class\n";
                    }
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    KRATOS_ERROR_IF(i_prototype == msRegisteredObjects.end())
                        << "There is no object registered in Kratos with name : "
                        << object_name << std::endl;

                    if(!pValue) {
                        pValue = std::move(Kratos::unique_ptr<TDataType>(static_cast<TDataType*>((i_prototype->second)())));
                    }
                }

                // Load the pointer address before loading the content
                mLoadedPointers[p_pointer]=pValue.get();
                load(rTag, *pValue);
            }
            else
            {
                pValue = std::move(Kratos::unique_ptr<TDataType>(static_cast<TDataType*>((i_pointer->second))));
            }
        }
    }

    template<class TDataType>
    void load(std::string const & rTag, TDataType*& pValue)
    {
        PointerType pointer_type = SP_INVALID_POINTER;
        void* p_pointer;
        read(pointer_type);

        if(pointer_type != SP_INVALID_POINTER)
        {
            read(p_pointer);
            LoadedPointersContainerType::iterator i_pointer = mLoadedPointers.find(p_pointer);
            if(i_pointer == mLoadedPointers.end())
            {
                if(pointer_type == SP_BASE_CLASS_POINTER)
                {
                    if constexpr (!std::is_abstract_v<TDataType>) {
                        if(!pValue) {
                            pValue = new TDataType;
                        }
                    }
                    else {
                        KRATOS_ERROR << "Cannot instantiate an abstract class\n";
                    }
                }
                else if(pointer_type == SP_DERIVED_CLASS_POINTER)
                {
                    std::string object_name;
                    read(object_name);
                    typename RegisteredObjectsContainerType::iterator i_prototype =  msRegisteredObjects.find(object_name);

                    KRATOS_ERROR_IF(i_prototype == msRegisteredObjects.end())
                        << "There is no object registered in Kratos with name : "
                        << object_name << std::endl;

                    if(!pValue) {
                        pValue = static_cast<TDataType*>((i_prototype->second)());
                    }

                }

                // Load the pointer address before loading the content
                mLoadedPointers[p_pointer]=&pValue;
                load(rTag, *pValue);
            }
            else
            {
                pValue = *static_cast<TDataType**>((i_pointer->second));
            }
        }
    }

    void load(std::string const & rTag, ModelPart*& pValue);

    void load(std::string const & rTag, Kratos::unique_ptr<ModelPart>& pValue);

    void load(std::string const & rTag, Kratos::shared_ptr<ModelPart>& pValue);


    template<class TDataType>
    void load(std::string const & rTag, Kratos::weak_ptr<TDataType>& pValue)
    {
    }

    template<class TDataType>
    void load(std::string const & rTag, WeakPointerVector<TDataType>& pValue)
    {
    }

    template<class TDataType>
    void load(std::string const & rTag, const Variable<TDataType>* /*pVariable*/)
    {
        load_trace_point(rTag);
        std::string name;
        read(name);
    }

    template<class TDataType, std::size_t TDataSize>
    void load(std::string const & rTag, std::array<TDataType, TDataSize>& rObject)
    {
        load_trace_point(rTag);
        for (SizeType i = 0; i < TDataSize; i++)
            load("E", rObject[i]);
    }

    template<class TDataType>
    void load(std::string const & rTag, std::vector<TDataType>& rObject)
    {
        load_trace_point(rTag);
        SizeType size;

        load("size", size);

        rObject.resize(size);

        for(SizeType i = 0 ; i < size ; i++)
            load("E", rObject[i]);
    }

    template<class TDataType>
    void load(std::string const & rTag, DenseVector<TDataType>& rObject)
    {
        load_trace_point(rTag);
        SizeType size;

        load("size", size);

        rObject.resize(size,false);

        for(SizeType i = 0 ; i < size ; i++)
            load("E", rObject[i]);
    }


    template<class TKeyType, class TDataType>
    void load(std::string const & rTag, std::map<TKeyType, TDataType>& rObject)
    {
        load_associative_container(rTag, rObject);
    }

    template<class TKeyType, class TDataType>
    void load(std::string const & rTag, std::unordered_map<TKeyType, TDataType>& rObject)
    {
        load_associative_container(rTag, rObject);
    }

    template<class TDataType>
    void load(std::string const & rTag, std::set<TDataType>& rObject)
    {
        load_associative_container(rTag, rObject);
    }

    template<class TDataType>
    void load(std::string const & rTag, std::unordered_set<TDataType>& rObject)
    {
        load_associative_container(rTag, rObject);
    }

    template<class TDataType, std::size_t TDimension>
    void load(std::string const & rTag, array_1d<TDataType, TDimension>& rObject)
    {
        load_trace_point(rTag);
        for(SizeType i = 0 ; i < TDimension ; i++)
            load("E", rObject[i]);
    }

    template<class TFirstType, class TSecondType>
    void load(std::string const & rTag, std::pair<TFirstType, TSecondType>& rObject)
    {
        load_trace_point(rTag);
        load("First", rObject.first);
        load("Second", rObject.second);
    }

    template<class TDataType, std::size_t TDimension>
    void load(std::string const & rTag, BoundedVector<TDataType, TDimension>& rObject)
    {
        load_trace_point(rTag);

        for(SizeType i = 0 ; i < TDimension ; ++i)
            load("E", rObject[i]);
    }

    template<class TDataType, std::size_t TDimension1, std::size_t TDimension2>
    void load(std::string const & rTag, BoundedMatrix<TDataType, TDimension1, TDimension2>& rObject)
    {
        load_trace_point(rTag);

        for(SizeType i = 0 ; i < TDimension1 ; ++i)
            for(SizeType j = 0 ; j < TDimension2 ; ++j)
                load("E", rObject(i,j));
//    read(rObject);
    }

    KRATOS_SERIALIZATION_DIRECT_LOAD(bool)
    KRATOS_SERIALIZATION_DIRECT_LOAD(int)
    KRATOS_SERIALIZATION_DIRECT_LOAD(long)
    KRATOS_SERIALIZATION_DIRECT_LOAD(double)
    KRATOS_SERIALIZATION_DIRECT_LOAD(float)
    KRATOS_SERIALIZATION_DIRECT_LOAD(unsigned long)
    KRATOS_SERIALIZATION_DIRECT_LOAD(unsigned int)
    KRATOS_SERIALIZATION_DIRECT_LOAD(std::string)
    KRATOS_SERIALIZATION_DIRECT_LOAD(Matrix)
    KRATOS_SERIALIZATION_DIRECT_LOAD(long long)
//#ifdef  _WIN32 // work around for windows int64_t error
//    KRATOS_SERIALIZATION_DIRECT_LOAD(__int64)
//#endif
#ifdef  _WIN64 // work around for windows size_t error in win64
    KRATOS_SERIALIZATION_DIRECT_LOAD(std::size_t)
#endif
    KRATOS_SERIALIZATION_DIRECT_LOAD(std::complex<double>)
    KRATOS_SERIALIZATION_DIRECT_LOAD(std::complex<float>)

    template<class TDataType, std::size_t TDataSize>
    void save(std::string const & rTag, std::array<TDataType, TDataSize> const& rObject)
    {
        save_trace_point(rTag);
        for (SizeType i = 0; i < TDataSize; i++)
            save("E", rObject[i]);
    }

    template<class TDataType>
    void save(std::string const & rTag, std::vector<TDataType> const& rObject)
    {
        save_trace_point(rTag);
        SizeType size = rObject.size();

        save("size", size);

        // Workaround for the special case of std::vector<bool>.
        // Long story short, the return type of std::vector<bool>::operator[]
        // is not necessarily bool&, but up to the STL vendor's decision.
        // Detailed explanation: https://github.com/KratosMultiphysics/Kratos/issues/10357#issuecomment-1274725614.
        using SaveType = std::conditional_t<
            std::is_same_v<typename std::decay<TDataType>::type, bool>,
            bool,
            const TDataType&
        >;

        for(SizeType i = 0 ; i < size ; i++)
            save("E", SaveType(rObject[i]));
    }

    template<class TDataType>
    void save(std::string const & rTag, DenseVector<TDataType> const& rObject)
    {
        save_trace_point(rTag);
        SizeType size = rObject.size();

        save("size", size);

        for(SizeType i = 0 ; i < size ; i++)
            save("E", rObject[i]);
    }

    template<class TDataType, std::size_t TDimension>
    void save(std::string const & rTag, array_1d<TDataType, TDimension> const& rObject)
    {
        save_trace_point(rTag);
        for(SizeType i = 0 ; i < TDimension ; i++)
            save("E", rObject[i]);
    }

    template<class TKeyType, class TDataType>
    void save(std::string const & rTag, std::map<TKeyType, TDataType> const& rObject)
    {
        save_associative_container(rTag, rObject);
    }

    template<class TKeyType, class TDataType>
    void save(std::string const & rTag, std::unordered_map<TKeyType, TDataType> const& rObject)
    {
        save_associative_container(rTag, rObject);
    }

    template<class TDataType>
    void save(std::string const & rTag, std::set<TDataType> const& rObject)
    {
        save_associative_container(rTag, rObject);
    }

    template<class TDataType>
    void save(std::string const & rTag, std::unordered_set<TDataType> const& rObject)
    {
        save_associative_container(rTag, rObject);
    }

    template<class TDataType>
    void save(std::string const & rTag, TDataType const& rObject)
    {
        save_trace_point(rTag);
        rObject.save(*this);
    }

    template<class TDataType>
    void save(std::string const & rTag, const Variable<TDataType>* pVariable)
    {
        save_trace_point(rTag);
        write(pVariable->Name());
    }


    template<class TDataType>
    void save(std::string const & rTag, Kratos::shared_ptr<TDataType> pValue)
    {
        save(rTag, pValue.get());
    }

    template<class TDataType>
    void save(std::string const & rTag, Kratos::intrusive_ptr<TDataType> pValue)
    {
        save(rTag, pValue.get());
    }


    template<class TDataType>
    void save(std::string const & rTag, Kratos::unique_ptr<TDataType> const& pValue)
    {
        save(rTag, pValue.get());
    }

    template<class TDataType>
    void save(std::string const & rTag, const TDataType * pValue)
    {
        if(pValue)
        {
            if(IsDerived(pValue))
                write(SP_DERIVED_CLASS_POINTER);
            else
                write(SP_BASE_CLASS_POINTER);

            SavePointer(rTag,pValue);
        }
        else
        {
            write(SP_INVALID_POINTER);
        }
    }

    template<class TDataType>
    bool IsDerived(TDataType * pValue)
    {
      if (strcmp(typeid(TDataType).name(), typeid(*pValue).name()) != 0) {
        return true;
      }
      else {
        return false;
      }
    }


    template<class TDataType>
    void save(std::string const & rTag, TDataType * pValue)
    {
        if(pValue)
        {
            if(IsDerived(pValue))
            {
                write(SP_DERIVED_CLASS_POINTER);
            }
            else
            {
                write(SP_BASE_CLASS_POINTER);
            }

            SavePointer(rTag,pValue);
        }
        else
        {
            write(SP_INVALID_POINTER);
        }
    }

    template<class TDataType>
    void save(std::string const & rTag, Kratos::weak_ptr<TDataType> pValue)
    {
    }

    template<class TDataType>
    void save(std::string const & rTag, Kratos::WeakPointerVector<TDataType> pValue)
    {
    }

    template<class TDataType>
    void save(std::string const & rTag, Kratos::shared_ptr<const TDataType> pValue)
    {
        save(rTag, pValue.get());
    }

    void save(std::string const & rTag, const char * pValue)
    {
        save_trace_point(rTag);
        write(std::string(pValue));
    }


    template<class TFirstType, class TSecondType>
    void save(std::string const & rTag, std::pair<TFirstType, TSecondType> rObject)
    {
        save_trace_point(rTag);
        save("First", rObject.first);
        save("Second", rObject.second);
    }

    template<class TDataType, std::size_t TDimension>
    void save(std::string const & rTag, BoundedVector<TDataType, TDimension> const& rObject)
    {
        save_trace_point(rTag);

        for(SizeType i = 0 ; i < TDimension ; ++i)
            save("E", rObject[i]);
    }

    template<class TDataType, std::size_t TDimension1, std::size_t TDimension2>
    void save(std::string const & rTag, BoundedMatrix<TDataType, TDimension1, TDimension2> const& rObject)
    {
        save_trace_point(rTag);

        for(SizeType i = 0 ; i < TDimension1 ; ++i)
            for(SizeType j = 0 ; j < TDimension2 ; ++j)
                save("E", rObject(i,j));
    }

    KRATOS_SERIALIZATION_DIRECT_SAVE(bool)
    KRATOS_SERIALIZATION_DIRECT_SAVE(int)
    KRATOS_SERIALIZATION_DIRECT_SAVE(long)
    KRATOS_SERIALIZATION_DIRECT_SAVE(double)
    KRATOS_SERIALIZATION_DIRECT_SAVE(float)
    KRATOS_SERIALIZATION_DIRECT_SAVE(unsigned long)
    KRATOS_SERIALIZATION_DIRECT_SAVE(unsigned int)
    KRATOS_SERIALIZATION_DIRECT_SAVE(std::string)
    KRATOS_SERIALIZATION_DIRECT_SAVE(Matrix)
    KRATOS_SERIALIZATION_DIRECT_SAVE(long long)
//#ifdef  _WIN32 // work around for windows int64_t error
//    KRATOS_SERIALIZATION_DIRECT_SAVE(__int64)
//#endif
#ifdef  _WIN64 // work around for windows size_t error in win64
    KRATOS_SERIALIZATION_DIRECT_SAVE(std::size_t)
#endif
    KRATOS_SERIALIZATION_DIRECT_SAVE(std::complex<double>)
    KRATOS_SERIALIZATION_DIRECT_SAVE(std::complex<float>)


    template<class TDataType>
    void load_base(std::string const & rTag, TDataType& rObject)
    {
        load_trace_point(rTag);
        rObject.TDataType::load(*this);
    }


    template<class TDataType>
    void load_base(std::string const & rTag, std::vector<TDataType>& rObject)
    {
        load_trace_point(rTag);
        load(rTag, rObject);
    }

    template<class TDataType>
    void load_base(std::string const & rTag, DenseVector<TDataType>& rObject)
    {
        load_trace_point(rTag);
        load(rTag, rObject);
    }

    template<class TDataType, std::size_t TDimension>
    void load_base(std::string const & rTag, array_1d<TDataType, TDimension>& rObject)
    {
        load_trace_point(rTag);
        load(rTag, rObject);
    }

    template<class TDataType>
    void save_base(std::string const & rTag, std::vector<TDataType> const& rObject)
    {
        save_trace_point(rTag);
        save(rTag, rObject);
    }

    template<class TDataType>
    void save_base(std::string const & rTag, DenseVector<TDataType> const& rObject)
    {
        save_trace_point(rTag);
        save(rTag, rObject);
    }

    template<class TDataType, std::size_t TDimension>
    void save_base(std::string const & rTag, array_1d<TDataType, TDimension> const& rObject)
    {
        save_trace_point(rTag);
        save(rTag, rObject);
    }

    template<class TDataType>
    void save_base(std::string const & rTag, TDataType const& rObject)
    {
        save_trace_point(rTag);
        rObject.TDataType::save(*this);
    }

    void save_trace_point(std::string const & rTag)
    {
        if(mTrace)
        {
            write(rTag);
        }
    }

    bool load_trace_point(std::string const & rTag)
    {
        if(mTrace == SERIALIZER_TRACE_ERROR) // only reporting the errors
        {
            std::string read_tag;
            read(read_tag);
            if(read_tag == rTag)
                return true;
            else
            {
                std::stringstream buffer;
                buffer << "In line " << mNumberOfLines;
                buffer << " the trace tag is not the expected one:" << std::endl;
                buffer << "    Tag found : " << read_tag << std::endl;
                buffer << "    Tag given : " << rTag << std::endl;
                KRATOS_ERROR << buffer.str() << std::endl;
            }
        }
        else if(mTrace == SERIALIZER_TRACE_ALL) // also reporting matched tags.
        {
            std::string read_tag;
            read(read_tag);
            if(read_tag == rTag)
            {
                KRATOS_INFO("Serializer") << "In line " << mNumberOfLines << " loading " << rTag << " as expected" << std::endl;
                return true;
            }
            else
            {
                std::stringstream buffer;
                buffer << "In line " << mNumberOfLines;
                buffer << " the trace tag is not the expected one:" << std::endl;
                buffer << "    Tag found : " << read_tag << std::endl;
                buffer << "    Tag given : " << rTag << std::endl;
                KRATOS_ERROR << buffer.str() << std::endl;
            }
        }
        return false;

    }




    ///@}
    ///@name Access
    ///@{

    BufferType* pGetBuffer()
    {
        return mpBuffer;
    }

    /**
     * This function let's one introduce "pValue"  between the objects
     * which are considered to be already serialized
     * TODO: verify if this should be a void* or if it is correct that it is taken as TDataType
     */
    template<class TDataType>
    void AddToSavedPointers(const TDataType& pValue) {
        mSavedPointers.insert(pValue);
    }

    /**
     * This function is to be used to inform the serializer that the object
     * initially stored in "pStoredPosition" is after loading located at pAllocatedPosition
     * This function is useful to substitute some objects with others that already exist.
     */
    void RedirectLoadingPointer(void * pStoredPointer, void * pAllocatedPosition) {
        mLoadedPointers[pStoredPointer]=pAllocatedPosition;
    }

    static RegisteredObjectsContainerType& GetRegisteredObjects()
    {
        return msRegisteredObjects;
    }

    static RegisteredObjectsNameContainerType& GetRegisteredObjectsName()
    {
        return msRegisteredObjectsName;
    }

    void Set(const Flags ThisFlag)
    {
        mFlags.Set(ThisFlag);
    }

    ///@}
    ///@name Inquiry
    ///@{

    bool Is(Flags const & rOtherFlag) const
    {
        return mFlags.Is(rOtherFlag);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Serializer";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static RegisteredObjectsContainerType msRegisteredObjects;
    static RegisteredObjectsNameContainerType msRegisteredObjectsName;

    ///@}
    ///@name Member Variables
    ///@{

    Flags mFlags;

    BufferType* mpBuffer;
    TraceType mTrace;
    SizeType mNumberOfLines;

    SavedPointersContainerType mSavedPointers;
    LoadedPointersContainerType mLoadedPointers;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    template<class TDataType>
    void SavePointer(std::string const & rTag, const TDataType * pValue)
    {
        write(pValue);
        if (mSavedPointers.find(pValue) == mSavedPointers.end())
        {
            mSavedPointers.insert(pValue);
            if (IsDerived(pValue))
            {
                typename RegisteredObjectsNameContainerType::iterator i_name = msRegisteredObjectsName.find(typeid (*pValue).name());

                if (i_name == msRegisteredObjectsName.end()) {
                    KRATOS_ERROR << "There is no object registered in Kratos with type id : "
                                 << typeid (*pValue).name() << std::endl;
                } else {
                    write(i_name->second);
                }
            }

            save(rTag, *pValue);
        }
    }

    VariableData* GetVariableData(std::string const & VariableName);

    template<class TMapType>
    void load_associative_container(std::string const & rTag, TMapType& rObject)
    {
        load_trace_point(rTag);
        SizeType size = rObject.size();

        load("size", size);

        for(SizeType i = 0 ; i < size ; i++){
            typename TMapType::value_type temp;
            load("E", temp);
            rObject.insert(temp);
        }
    }

    template<class TMapType>
    void save_associative_container(std::string const & rTag, TMapType const& rObject)
    {
        save_trace_point(rTag);
        SizeType size = rObject.size();

        save("size", size);

        for(auto& i : rObject)
            save("E", i);
    }

    void read(PointerType& rValue)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        int temp;
        mpBuffer->read((char *)(&temp),sizeof(PointerType));
        rValue = PointerType(temp);

        KRATOS_SERIALIZER_MODE_ASCII

        int temp;
        *mpBuffer >> temp;
        rValue = PointerType(temp);
        mNumberOfLines++;

        KRATOS_SERIALIZER_MODE_END
    }

    void write(PointerType const& rValue)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        int ptr = (int)rValue;
        const char * data = reinterpret_cast<const char*>(&ptr);
        mpBuffer->write(data,sizeof(PointerType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << int(rValue) << std::endl;

        KRATOS_SERIALIZER_MODE_END
    }

    void read(std::string& rValue)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType size;
        mpBuffer->read((char *)(&size),sizeof(SizeType));
        rValue.resize(size);
        if (size>0) {
            mpBuffer->read(rValue.data(), size);
        }

        KRATOS_SERIALIZER_MODE_ASCII

        // going to the first '"'
        std::getline( *mpBuffer,rValue, '\"');
        // reading the string itself until second '"'
        std::getline( *mpBuffer,rValue, '\"');
        mNumberOfLines++;

        KRATOS_SERIALIZER_MODE_END
    }

    void write(std::string const& rValue)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        const char * data = rValue.c_str();
        SizeType rData_size = rValue.length() * sizeof(char);

        const char * data1 = reinterpret_cast<const char *>(&rData_size);

        mpBuffer->write(data1,sizeof(SizeType));
        mpBuffer->write(data,rData_size);

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << "\"" << rValue << "\"" << std::endl;

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void read(TDataType& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        mpBuffer->read((char *)(&rData),sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer >> rData;
        mNumberOfLines++;

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void write(TDataType const& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        const char * data = reinterpret_cast<const char*>(&rData);
        mpBuffer->write(data,sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << rData << std::endl;

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void read(std::vector<TDataType>& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType size;
        mpBuffer->read((char *)(&size),sizeof(SizeType));

        rData.resize(size);

        read(rData.begin(), rData.end(), sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        std::size_t size;
        *mpBuffer >> size;
        rData.resize(size);
        mNumberOfLines++;

        read(rData.begin(), rData.end());

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void write(std::vector<TDataType> const& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType rData_size = rData.size();

        const char * data = reinterpret_cast<const char *>(&rData_size);
        mpBuffer->write(data,sizeof(SizeType));

        write(rData.begin(), rData.end(), sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << rData.size() << std::endl;
        write(rData.begin(), rData.end());

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void read(DenseMatrix<TDataType>& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType size1;
        SizeType size2;

        mpBuffer->read((char *)(&size1),sizeof(SizeType));
        mpBuffer->read((char *)(&size2),sizeof(SizeType));

        rData.resize(size1,size2);

        read(rData.data().begin(), rData.data().end(), sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        SizeType size1;
        SizeType size2;

        *mpBuffer >> size1;
        mNumberOfLines++;
        *mpBuffer >> size2;
        mNumberOfLines++;

        rData.resize(size1,size2);

        read(rData.data().begin(), rData.data().end(),0);

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TDataType>
    void write(DenseMatrix<TDataType> const& rData)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        SizeType rData_size1 = rData.size1();
        SizeType rData_size2 = rData.size2();

        const char * data1 = reinterpret_cast<const char *>(&rData_size1);
        const char * data2 = reinterpret_cast<const char *>(&rData_size2);

        mpBuffer->write(data1,sizeof(SizeType));
        mpBuffer->write(data2,sizeof(SizeType));

        write(rData.data().begin(), rData.data().end(), sizeof(TDataType));

        KRATOS_SERIALIZER_MODE_ASCII

        *mpBuffer << rData.size1() << std::endl;
        *mpBuffer << rData.size2() << std::endl;

        write(rData.data().begin(), rData.data().end(),0);

        KRATOS_SERIALIZER_MODE_END
    }

    template<class TIteratorType>
    void read(TIteratorType First, TIteratorType Last, SizeType size)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        for(; First != Last ; First++)
        {
            mpBuffer->read((char *)First,sizeof(size));
        }

        KRATOS_SERIALIZER_MODE_ASCII

        for(; First != Last ; First++)
        {
            *mpBuffer >> *First;
            mNumberOfLines++;

        }

        KRATOS_SERIALIZER_MODE_END
    }
    template<class TIteratorType>
    void write(TIteratorType First, TIteratorType Last, SizeType size)
    {
        KRATOS_SERIALIZER_MODE_BINARY

        for(; First != Last ; First++)
        {
            const char * data = reinterpret_cast<const char *>(First);
            mpBuffer->write(data,sizeof(size));
        }

        KRATOS_SERIALIZER_MODE_ASCII

        for(; First != Last ; First++)
            *mpBuffer << *First << std::endl;

        KRATOS_SERIALIZER_MODE_END
    }

    inline SizeType BlockCompatibleSize(SizeType rSize)
    {
        typedef char BlockType;
        const SizeType block_size = sizeof(BlockType);
        return static_cast<SizeType>(((block_size - 1) + rSize) / block_size);
    }

    /// Sets the pointer of the stream buffer at the begnining
    void SeekBegin();

    /// Sets the pointer of the stream buffer at the end
    void SeekEnd();

    template <typename T>
    std::shared_ptr<T> MakeShared()
    {
        if constexpr (std::is_default_constructible_v<T>)
        {
            return std::make_shared<T>();
        } else {
            return std::shared_ptr<T>(new T);
        }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Serializer& operator=(Serializer const& rOther);

    /// Copy constructor.
    Serializer(Serializer const& rOther);

    ///@}

}; // Class Serializer

///@}

}  // namespace Kratos.

#undef KRATOS_SERIALIZER_MODE_BINARY
#undef KRATOS_SERIALIZER_MODE_ASCII
#undef KRATOS_SERIALIZER_MODE_END

#undef KRATOS_SERIALIZATION_DIRECT_LOAD
#undef KRATOS_SERIALIZATION_DIRECT_SAVE
