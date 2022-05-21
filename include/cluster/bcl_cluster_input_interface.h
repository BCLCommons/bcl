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

#ifndef BCL_CLUSTER_INPUT_INTERFACE_H_
#define BCL_CLUSTER_INPUT_INTERFACE_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_ifstream.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_object_nd_hash_map.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class InputInterface
    //! @brief InputInterface is the interface class from which methods for inputting from different file formats
    //!        should be derived.
    //!
    //! @remarks example unnecessary
    //! @author alexanns, woetzen
    //! @date August 20, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class InputInterface :
      public util::ObjectInterface
    {
    protected:

    //////////
    // data //
    //////////

      //! m_Objects is a ShPtr to the list of actual objects to be clustered
      util::ShPtr< storage::List< t_DataType> > m_Objects;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      InputInterface() :
        m_Objects()
      {
      }

      //! @brief constructor from the objects
      //! @param OBJECTS is string denoting the file containing the data which will be "m_Objects"
      InputInterface
      (
        const util::ShPtr< storage::List< t_DataType> > &OBJECTS
      ) :
        m_Objects( OBJECTS)
      {
      }

      //! @brief copy constructor
      //! @param INPUT_INTERFACE the InputInterface from which this InputInterface will be copied
      InputInterface( const InputInterface &INPUT_INTERFACE) :
        m_Objects( INPUT_INTERFACE.m_Objects)
      {
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief GetInputObjects gives the objects which have been created
      //! @return util::ShPtr< storage::List< t_DataType> > which are the objects that were created
      const util::ShPtr< storage::List< t_DataType> > &GetInputObjects() const
      {
        return m_Objects;
      }

      //! @brief SetInputObjects gives the objects which have been created
      //! param OBJECTS the objects that m_Objects will be set to
      void SetInputObjects
      (
        const util::ShPtr< storage::List< t_DataType> > &OBJECTS
      )
      {
        m_Objects = OBJECTS;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief HandleInput gets data from input; puts it into the data construct with pointers to the objects
      //         in m_Objects
      //! @param IFSTREAM is the stream from which the input will be read
      //! @return returns data construct holding the distances between all objects for use in a LinkageInterface
      virtual
      util::ShPtr
      <
        math::FunctionInterface
        <
          storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType
        >
      >
      HandleInput( io::IFStream &IFSTREAM) = 0;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_Objects, ISTREAM);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Objects, OSTREAM, INDENT);

        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // InputInterface

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_INPUT_INTERFACE_H_ 
