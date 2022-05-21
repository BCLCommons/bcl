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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_molecule_shape_moments.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeShapeMoments::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeShapeMoments()
      )
    );

    //! @brief default constructor
    MoleculeShapeMoments::MoleculeShapeMoments()
    {
    }

    //! @brief constructor from atom properties
    MoleculeShapeMoments::MoleculeShapeMoments
    (
      const CheminfoProperty &ATOM_PROPERTY_CENTERING,
      const CheminfoProperty &ATOM_PROPERTY_WEIGHTING
    ) :
      m_AtomPropertyAnchorPoints( ATOM_PROPERTY_CENTERING),
      m_AtomPropertyWeighting( ATOM_PROPERTY_WEIGHTING)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeShapeMoments
    MoleculeShapeMoments *MoleculeShapeMoments::Clone() const
    {
      return new MoleculeShapeMoments( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeShapeMoments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeShapeMoments::GetAlias() const
    {
      static const std::string s_name( "ShapeMoments");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > MoleculeShapeMoments::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomPropertyAnchorPoints,
          &m_AtomPropertyWeighting + 1
        );
    }

    namespace
    {
      //! @brief compute skew, and return weighted statistics
      float SkewAndStats
      (
        const linal::VectorConstInterface< float> &VEC,
        const linal::VectorConstInterface< float> &WEIGHTS,
        math::RunningAverageSD< float> &STATS_COLLECTOR
      )
      {
        STATS_COLLECTOR.Reset();
        for( size_t i( 0), sz( VEC.GetSize()); i < sz; ++i)
        {
          STATS_COLLECTOR.AddWeightedObservation( VEC( i), WEIGHTS( i));
        }
        math::RunningAverage< float> skew_ave;
        for( size_t i( 0), sz( VEC.GetSize()); i < sz; ++i)
        {
          skew_ave.AddWeightedObservation( math::Pow( VEC( i) - STATS_COLLECTOR.GetAverage(), float( 3.0)), WEIGHTS( i));
        }
        const float skew_std
        (
          std::max( STATS_COLLECTOR.GetSampleStandardDeviation(), std::numeric_limits< float>::epsilon())
        );
        const float skew
        (
          skew_ave.GetAverage() / math::Pow( skew_std, float( 3.0))
        );
        const float skew_weight( VEC.GetSize());
        return skew * skew_weight * skew_weight / ( std::max( skew_weight, float( 1.5)) - 1.0) / ( std::max( skew_weight, float( 2.5)) - 2.0);
      }

      //! @brief compute distances from a point on a molecule
      linal::Vector< float> GetDistancesTo
      (
        const linal::Vector3D &POINT,
        const SequenceInterface< chemistry::AtomConformationalInterface> &MOL
      )
      {
        linal::Vector< float> distances( MOL.GetSize());
        size_t pos( 0);
        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr( MOL.GetIterator());
          itr.NotAtEnd();
          ++itr, ++pos
        )
        {
          distances( pos) = linal::Distance( POINT, itr->GetPosition());
        }
        return distances;
      }
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void MoleculeShapeMoments::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      util::SiPtr< const coord::OrientationInterface> orientation
      (
        this->GetCurrentObject()
      );

      math::RunningAverage< linal::Vector3D> weighted_center;
      iterate::Generic< const chemistry::AtomConformationalInterface> itr_atom( this->GetCurrentObject()->GetIterator());
      linal::Vector< float> anchor_weights( itr_atom.GetSize()), distance_weights( itr_atom.GetSize());
      size_t position( 0);
      for
      (
        Iterator< chemistry::AtomConformationalInterface> itr_desc( itr_atom);
        itr_atom.NotAtEnd();
        ++itr_atom, ++itr_desc, ++position
      )
      {
        const float weight( anchor_weights( position) = ( *m_AtomPropertyAnchorPoints)( itr_desc)( 0));
        weighted_center.AddWeightedObservation( itr_atom->GetPosition(), weight);
        distance_weights( position) = ( *m_AtomPropertyWeighting)( itr_desc)( 0);
      }
      if( anchor_weights.Sum() <= float( 0.0) || distance_weights.Sum() <= float( 0.0))
      {
        STORAGE = float( 0.0);
        return;
      }
      linal::Vector< float> dist_to_center( GetDistancesTo( weighted_center.GetAverage(), *this->GetCurrentObject()));
      math::RunningAverageSD< float> dist_to_center_stats;
      float centroid_skew( SkewAndStats( dist_to_center, distance_weights, dist_to_center_stats));
      STORAGE( 0) = dist_to_center_stats.GetAverage();
      STORAGE( 1) = dist_to_center_stats.GetSampleStandardDeviation();
      STORAGE( 2) = centroid_skew;

      // closest index to centroid
      size_t closest_index( math::Statistics::MinimumIndex( dist_to_center.Begin(), dist_to_center.End()));
      size_t furthest_index( math::Statistics::MaximumIndex( dist_to_center.Begin(), dist_to_center.End()));
      itr_atom.GotoPosition( closest_index);
      dist_to_center = GetDistancesTo( itr_atom->GetPosition(), *this->GetCurrentObject());
      centroid_skew = SkewAndStats( dist_to_center, distance_weights, dist_to_center_stats);
      STORAGE( 3) = dist_to_center_stats.GetAverage();
      STORAGE( 4) = dist_to_center_stats.GetSampleStandardDeviation();
      STORAGE( 5) = centroid_skew;
      itr_atom.GotoPosition( furthest_index);
      dist_to_center = GetDistancesTo( itr_atom->GetPosition(), *this->GetCurrentObject());
      centroid_skew = SkewAndStats( dist_to_center, distance_weights, dist_to_center_stats);
      STORAGE( 6) = dist_to_center_stats.GetAverage();
      STORAGE( 7) = dist_to_center_stats.GetSampleStandardDeviation();
      STORAGE( 8) = centroid_skew;
      furthest_index = math::Statistics::MaximumIndex( dist_to_center.Begin(), dist_to_center.End());
      itr_atom.GotoPosition( furthest_index);
      dist_to_center = GetDistancesTo( itr_atom->GetPosition(), *this->GetCurrentObject());
      centroid_skew = SkewAndStats( dist_to_center, distance_weights, dist_to_center_stats);
      STORAGE( 9) = dist_to_center_stats.GetAverage();
      STORAGE( 10) = dist_to_center_stats.GetSampleStandardDeviation();
      STORAGE( 11) = centroid_skew;
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeShapeMoments::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes the mean, std, and skew of atomic distances from four anchor points: the molecular centroid (mct), "
        "closest atom to the mct (cam), furthest atom from mct (fam), and furthest atom from fam (faf). "
        "The molecular centroid is the center of the positions, weighted by the anchor property."
      );

      parameters.AddInitializer
      (
        "anchor",
        "weighting property used to define the centroid and other anchor points in the target molecule",
        io::Serialization::GetAgent( &m_AtomPropertyAnchorPoints)
      );
      parameters.AddInitializer
      (
        "weighting",
        "weighting for statistical moments; e.g. use Atom_Polarizability to place more weight on distances associated with "
        "highly polar species in the molecule",
        io::Serialization::GetAgent( &m_AtomPropertyWeighting)
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeShapeMoments::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // catch m_NumberSteps == 0
      if( m_AtomPropertyWeighting.IsDefined() && m_AtomPropertyWeighting->GetNormalSizeOfFeatures() != 1)
      {
        ERR_STREAM
           << "Expected a property that returned 1 properties per atom, but property returns "
           << m_AtomPropertyWeighting->GetNormalSizeOfFeatures()
           << " values per atom ( property was "
           << m_AtomPropertyWeighting->GetString()
           << ")";
        return false;
      }
      if( m_AtomPropertyAnchorPoints.IsDefined() && m_AtomPropertyAnchorPoints->GetNormalSizeOfFeatures() != 1)
      {
        ERR_STREAM
           << "Expected a property that returned 1 properties per atom, but property returns "
           << m_AtomPropertyAnchorPoints->GetNormalSizeOfFeatures()
           << " values per atom ( property was "
           << m_AtomPropertyAnchorPoints->GetString()
           << ")";
        return false;
      }
      m_AtomPropertyAnchorPoints->SetDimension( 1);
      m_AtomPropertyWeighting->SetDimension( 1);
      return true;
    }

  } // namespace descriptor
} // namespace bcl
