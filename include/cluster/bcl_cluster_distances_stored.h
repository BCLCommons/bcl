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

#ifndef BCL_CLUSTER_DISTANCES_STORED_H_
#define BCL_CLUSTER_DISTANCES_STORED_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_object_nd_hash_map.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DistancesStored
    //! @brief DistancesStored is a method for getting the distance between two objects which are being clustered.
    //! @details The    LinkageInterface asks a FunctionInterface for the distance between two objects when calculating the
    //!        linkage between two clusters. DistancesStored is one FunctionInterface which can be used to provide the
    //!        distance between two of the objects being clustered. It stores the distances in memory.
    //!
    //! @see @link example_cluster_distances_stored.cpp @endlink
    //! @author alexanns
    //! @date August 17, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class DistancesStored :
      public math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType>
    {
    private:

    //////////
    // data //
    //////////

      //! m_Data holds SiPtrs to objects and the distances between them
      storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > m_Data;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DistancesStored< t_DataType, t_PrecisionType>() :
        m_Data()
      {
      }

      //! @brief constructor from an object of the same type as "m_Data"
      DistancesStored< t_DataType, t_PrecisionType>
      (
        const storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > &DATA
      ) :
        m_Data( DATA)
      {
        // make sure all distances are defined
        for
        (
          typename storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> >::const_iterator
            itr_a( m_Data.Begin()), itr_a_end( m_Data.End());
          itr_a != itr_a_end;
          ++itr_a
        )
        {
          for
          (
            typename storage::HashMap< size_t, t_PrecisionType>::const_iterator
              itr_b( itr_a->second.Begin()), itr_b_end( itr_a->second.End());
            itr_b != itr_b_end;
            ++itr_b
          )
          {
            BCL_Assert( util::IsDefined( itr_b->second), "pairwise distance is undefined. are there nans in your input file?");
          }
        }
      }

      //! @brief copy constructor
      //! @param DISTANCE_STORED the DistancesStored< t_DataType> from which this DistancesStored will be copied
      DistancesStored< t_DataType, t_PrecisionType>( const DistancesStored< t_DataType, t_PrecisionType> &DISTANCE_STORED) :
        m_Data( DISTANCE_STORED.m_Data)
      {
      }

      //! @brief Clone function
      //! @return pointer to new DistancesStored< t_DataType>
      virtual DistancesStored< t_DataType, t_PrecisionType> *Clone() const
      {
        return new DistancesStored< t_DataType, t_PrecisionType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief GetData returns a const reference to "m_Data"
      //! @return returns m_Data const reference
      const storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > &
      GetData() const
      {
        return m_Data;
      }

      //! @brief SetData set "m_Data" to different data
      //! @param DATA is the new data that "m_Data" will be set to
      void SetData
      (
        const storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> > &DATA
      )
      {
        m_Data = DATA;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() which takes two t_DataTypes and returns the distance between them
      //! @param OBJECTS is a VectorND which has the two objects whose distance is needed
      //! @return returns a t_PrecisionType which is the distance between the two objects
      virtual t_PrecisionType operator()( const storage::VectorND< 2, util::SiPtr< const t_DataType> > &OBJECTS) const
      {
        //// create t_PrecisionType "distance" to hold the distance between the two objects of "OBJECTS"
        //t_PrecisionType distance;
        typename storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> >::const_iterator
          itr_find( m_Data.Find( size_t( OBJECTS.First().GetPointer())));

        const typename storage::HashMap< size_t, storage::HashMap< size_t, t_PrecisionType> >::const_iterator
          data_end_itr( m_Data.End());

        // this will only be true in the special case that first is the very last object in the object list
        // all other objects should be found able to be found as first key although the order of lookup may need to be
        // switched as handled further below
        if( itr_find == data_end_itr)
        {
          itr_find = m_Data.Find( size_t( OBJECTS.Second().GetPointer()));
          BCL_Assert( itr_find != data_end_itr, "could not find " + util::Format()( OBJECTS));
          const typename storage::HashMap< size_t, t_PrecisionType>::const_iterator itr_find_second( itr_find->second.Find( size_t( OBJECTS.First().GetPointer())));
          const typename storage::HashMap< size_t, t_PrecisionType>::const_iterator itr_find_second_end( itr_find->second.End());

          BCL_Assert( itr_find_second != itr_find_second_end, "could not find second in\n" + util::Format()( OBJECTS));

          return itr_find_second->second;
        }

        BCL_Assert( itr_find != data_end_itr, "could not find " + util::Format()( OBJECTS));

        typename storage::HashMap< size_t, t_PrecisionType>::const_iterator itr_find_second( itr_find->second.Find( size_t( OBJECTS.Second().GetPointer())));
        typename storage::HashMap< size_t, t_PrecisionType>::const_iterator itr_find_second_end( itr_find->second.End());

        // if true means that due to the upper triangle second first pairwise distance was calculated instead of
        // first second need to search for second first instead of first second
        if( itr_find_second == itr_find_second_end)
        {
          itr_find = m_Data.Find( size_t( OBJECTS.Second().GetPointer()));
          BCL_Assert( itr_find != data_end_itr, "could not find " + util::Format()( OBJECTS));
          itr_find_second = itr_find->second.Find( size_t( OBJECTS.First().GetPointer()));
          itr_find_second_end = itr_find->second.End();
        }

        BCL_Assert( itr_find_second != itr_find_second_end, "could not find second in\n" + util::Format()( OBJECTS));

        return itr_find_second->second;
      }

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
        io::Serialize::Read( m_Data, ISTREAM);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Data, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // DistancesStored

    // instantiate s_Instance
    template< typename t_DataType, typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> DistancesStored< t_DataType, t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new DistancesStored< t_DataType, t_PrecisionType>())
    );

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_DISTANCES_STORED_H_ 
