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

// include header of this class
#include "bcl_descriptor_atom_planarity.h"
// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_3d.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_linear_least_squares.h"
#include "math/bcl_math_running_average.h"
#include "type/bcl_type_compare.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief Clone function
    //! @return pointer to new AtomPlanarity
    template< typename t_DataType>
    AtomPlanarity< t_DataType> *AtomPlanarity< t_DataType>::Clone() const
    {
      return new AtomPlanarity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &AtomPlanarity< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &AtomPlanarity< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "PlanarityAtoms");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    //! @note only one value is needed if this is a numeric descriptor, for char descriptors, assume that 99999 is the
    //! @note max, so 5 characters is sufficient
    template< typename t_DataType>
    size_t AtomPlanarity< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer AtomPlanarity< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Returns the chi-squared value of the set of atoms from a perfect plane");
      parameters.AddInitializer
      (
        "",
        "coordinate retrievers; must each return 3 values per atom",
        io::Serialization::GetAgentWithSizeLimits( &m_CoordinateRetrievers, size_t( 4))
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool AtomPlanarity< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      for
      (
        typename storage::Vector< t_Implementation>::const_iterator
          itr( m_CoordinateRetrievers.Begin()), itr_end( m_CoordinateRetrievers.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetSizeOfFeatures() != size_t( 3))
        {
          ERR_STREAM
            << "AtomPlanarity requires descriptors that return coordinates; 3 values per atom, but "
            << ( *itr)->GetString() << " returns " << ( *itr)->GetSizeOfFeatures();
          return false;
        }
      }
      m_CoordinatesStorage = linal::Matrix< float>( m_CoordinateRetrievers.GetSize(), size_t( 3));
      m_CoordinatesStorage2D = linal::Matrix< float>( m_CoordinateRetrievers.GetSize(), size_t( 2));
      m_ZStorage = linal::Vector< float>( m_CoordinateRetrievers.GetSize(), float( 0.0));
      return true;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void AtomPlanarity< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // compute the centroid
      math::RunningAverage< linal::Vector< float> > centroid_ave;

      // create an iterator from this element
      Iterator< t_DataType> itr_element( ELEMENT);
      size_t row( 0);
      for
      (
        typename storage::Vector< t_Implementation>::iterator
          itr( m_CoordinateRetrievers.Begin()), itr_end( m_CoordinateRetrievers.End());
        itr != itr_end;
        ++itr, ++row
      )
      {
        // get coordinates from this
        const linal::VectorConstReference< float> coord( ( **itr)( itr_element));

        // add it to the centroid
        centroid_ave += coord;

        if( !coord.IsDefined())
        {
          STORAGE( 0) = util::GetUndefined< float>();
          return;
        }

        // insert it into the matrix
        std::copy( coord.Begin(), coord.End(), m_CoordinatesStorage[ row]);
      }

      // get the centroid x, y and z
      const linal::Vector< float> &centroid( centroid_ave.GetAverage());

      // subtract the centroid from all data points
      for( size_t row( 0), n_rows( m_CoordinateRetrievers.GetSize()); row < n_rows; ++row)
      {
        linal::VectorReference< float> rowref( size_t( 3), m_CoordinatesStorage[ row]);
        rowref -= centroid;
        rowref.Normalize();
      }
      BCL_MessageVrb( "Coordinates: " + util::Format()( m_CoordinatesStorage));
      linal::Vector< float> chixyz( size_t( 3));
      for( size_t dim( 0); dim < size_t( 3); ++dim)
      {
        size_t xdim1( dim == size_t( 0) ? 1 : 0);
        size_t xdim2( dim == size_t( 2) ? 1 : 2);
        for( size_t row( 0), n_rows( m_CoordinateRetrievers.GetSize()); row < n_rows; ++row)
        {
          m_CoordinatesStorage2D( row, 0) = m_CoordinatesStorage[ row][ xdim1];
          m_CoordinatesStorage2D( row, 1) = m_CoordinatesStorage[ row][ xdim2];
          m_ZStorage( row) = m_CoordinatesStorage[ row][ dim];
        }
        STORAGE( 0) += math::LinearLeastSquares::SolutionAndChiSquared( m_CoordinatesStorage2D, m_ZStorage).Second();
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > AtomPlanarity< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( m_CoordinateRetrievers.Begin(), m_CoordinateRetrievers.End());
    }

  } // namespace descriptor
} // namespace bcl
