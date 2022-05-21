// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_UTIL_ENUM_DATA_H_
#define BCL_UTIL_ENUM_DATA_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_class_descriptor.h"
#include "bcl_util_undefined.h"
#include "io/bcl_io_serialize.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnumData
    //! @brief gathers the index, name and the t_DataType in one object for later storage in sets
    //! @details EnumData is derived from t_DataType class and adds and Index and a Name
    //!
    //! @tparam t_DataType the class the is to be enumerated. t_DataType has a Read and Write function,
    //!         default and copy constructor
    //!
    //! @see @link example_util_enum_data.cpp @endlink
    //! @author heinzes1, woetzen, karakam
    //! @date Nov 4, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class EnumData :
      public t_DataType
    {

    public:

    //////////
    // data //
    //////////

      size_t      m_Index; //! index of enumerated object
      std::string m_Name;  //! name of enumerated object

      //! @brief returns name for undefined Enum
      //! @return reference to string that is Undefined Enum name
      static const std::string &GetUndefinedEnumName()
      {
        // static string to store undefined name
        static const std::string s_undefined_enum_name( "Undefined");

        // return
        return s_undefined_enum_name;
      }

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct undefined enum object
      EnumData() :
        t_DataType(),
        m_Index( GetUndefinedSize_t()),
        m_Name( GetUndefinedEnumName())
      {
      }

      //! @brief construct enum object - only used by this class, so that INDEX always has the proper numbering
      //! @param INDEX  index of current object
      //! @param NAME   name  of current object
      //! @param OBJECT          current enumerated object
      EnumData
      (
        const size_t INDEX,
        const std::string &NAME,
        const t_DataType &OBJECT
      ) :
        t_DataType( OBJECT),
        m_Index( INDEX),
        m_Name( NAME)
      {
      }

    ///////////////
    // operators //
    ///////////////

      //! get integer representation of enum data object
      //! @return index of current enum
      size_t GetIndex() const
      {
        return m_Index;
      }

      //! get associated object of enum data object
      //! @return name of current enum
      const std::string &GetName() const
      {
        return m_Name;
      }

      //! @brief assign enum data from given EnumDataObject
      //! @details does not assign the index or name
      //! @return reference to this
      EnumData< t_DataType> &operator =( const EnumData< t_DataType> &ENUM_DATA_RHS)
      {
        // just assign the t_DataType
        t_DataType::operator=( ENUM_DATA_RHS);

        // end
        return *this;
      }

      //! @brief singleton for undefined EnumData
      //! @details The undefined enum data does not belong to any particular derived Enum class.
      //! @return const reference to only instance of the undefined data
      static const EnumData< t_DataType> &GetUndefinedData()
      {
        static const EnumData s_undefined_data;

        return s_undefined_data;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write to std::ostream, abstract, depending of template parameter
      //! @param OSTREAM stream to write to
      //! @param INDENT indentation of output
      //! @return std::ostream the stream written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_Index, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Name, OSTREAM, INDENT) << '\n';

        // write base
        t_DataType::Write( OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

      //! @brief read from std::istream, abstract, depending of template parameter
      //! @param ISTREAM stream to read from
      //! @return std::istream the stream written to
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_Index, ISTREAM);
        io::Serialize::Read( m_Name, ISTREAM);

        // read base
        t_DataType::Read( ISTREAM);

        // end
        return ISTREAM;
      }

    }; // template class EnumData

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_ENUM_DATA_H_
