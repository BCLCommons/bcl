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
#include "quality/bcl_quality_gdt.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "quality/bcl_quality_average.h"
#include "quality/bcl_quality_lcs.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
//    http://predictioncenter.org/casp/casp7/public/doc/LCS_GDT.README
//        LCS and GDT description
//
//    Longest Continuous Segments under specified CA RMSD cutoff (LCS).
//    The algorithm identifies all the longest continuous segments of residues
//    in the prediction deviating from the target by not more than specified
//    CA RMSD cutoff using many different superpositions.
//    Each residue in a prediction is assigned to the longest of such segments
//    provided if is a part of that segment. The absolutely longest continuous
//    segment in prediction under given RMSD cutoff is reported as well.
//    For different values of the CA RMSD cutoff (1.0 A, 2.0 A, and 5.0 A) the results
//    of the analysis are reported.
//    This measure can be used to evaluate ab initio 3D and comparative modeling
//    predictions.
//
//    Global Distance Test (GDT). The algorithm identifies in the prediction
//    the sets of residues deviating from the target by not more than
//    specified CA DISTANCE cutoff using many different superpositions.
//    Each residue in a prediction is assigned to the largest set of the residues
//    (not necessary continuous) deviating from the target by no more than a
//    specified distance cutoff.
//    This measure can be used to evaluate ab-initio 3D and comparative modeling
//    predictions.
//    For different values of DISTANCE cutoff (0.5 A, 1.0 A, 1.5 A, ... 10.0 A),
//    several measures are reported:
//    NUMBER_OF_CA_max  - the number of CA's from the "largest set" that
//           can fit under specified distance cutoff
//    PERCENT_OF_CA_Tg  - percent of CA's from the "largest set" comparing
//           to the total number of CA's in target
//    FRAGMENT: Beg-End - beginning and end of the segment containing the
//           "largest set" of CA's
//    RMS_LOCAL         - RMSD (root mean square deviation) calculated on the
//           "largest set" of CA's
//    RMS_ALL_CA        - RMSD calculated on all CA after superposition of
//           the prediction structure to the target structure
//           based on the "largest set" of CA's
//
//    The goal of introducing these two measures (GDT and LCS) is to provide a
//    tool that can be used for better detection of relatively good or bad parts
//    of the model.
//    - Using LCS we can localize the "best" continuous (along
//    the sequence) parts of the model that can fit under
//    RMSD thresholds: 1A, 2A, and 5A
//    Three blue lines represent the longest continuous sets
//    of residues that can fit under 1A, 2A, and 5A cutoff,
//    respectively.
//    - Using GDT we can localize the "best" sets of residues
//    (not necessary continuous) that can fit under DISTANCE
//    thresholds: 0.5A, 1.0A, 1.5A ,..., 10.0A
//    There are three blue lines on the GDT plot.
//    Each line represents the set of 5, 10, or 50 percent of
//    residues that can fit under specific distance cutoff (axis Y).
//    So, the lowest line represents residues (axis X) from the 5
//    percent sets of all target residues. Middle line identifies
//    those residues from the 10 percent sets, and highest from
//    50 percent sets.
//
//    The differences between LCS and GDT are the following:
//
//    1) LCS (Longest Continuous Segment) is based on RMSD cutoff.
//    2) The goal of LCS is to localize the longest continuous segment
//    of residues that can fit under RMSD cutoff.
//    3) Each residue in a prediction is assigned to the longest continuous
//    segment provided if is a part of that segment.
//    4) The data provided in the result files contains the LCS calculated
//    under three selected values of CA RMSD cutoff: 1A, 2A, and 5A
//
//    5) GDT (Global Distance Test) is based on the DISTANCE cutoff.
//    6) The goal of GDT is to localize the largest set of residues
//    (not necessary continuous) deviating from the target by no more than
//    a specified DISTANCE cutoff.
//    7) Each residue in a prediction is assigned to the largest set of the
//    residues provided if is a part of that set.
//    8) The data provided in the result files contains the GDT calculated under
//    several values of DISTANCE cutoff: 0.5, 1.0, 1.5, ... , 10.0 Angstroms.
//
//    Results of the analysis given by LCS algorithm show rather local features of
//    the model, while the residues considered in GDT come from the whole model
//    structure (they do not have to maintain the continuity along the sequence).
//
//    The GDT procedure is the following. Each three-residue segment and each
//    continuous segment found by LCS is used as a starting point to give an
//    initial equivalencies (model-target CA pairs) for a superposition.
//    The list of equivalencies is iteratively extended to produce the largest
//    set of residues that can fit under considered distance cutoff.
//    For collecting data about largest sets of residues the iterative
//    superposition procedure (ISP) is used.
//    The goal of the ISP method is to exclude from the calculations atoms
//    that are more than some threshold (cutoff) distance between the
//    model and the target structure after the transform is applied.
//    Starting from the initial set of atoms (C-alphas) the algorithm is the
//    following:
//    a) obtain the transform
//    b) apply the transform
//    c) identify all atom pairs for which distance is larger than the
//    threshold
//    d) re-obtain the transform, excluding those atoms
//    e) repeat b) - d) until the set of atoms used in calculations
//    is the same for two cycles running
//
//    -------

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> GDT::s_Instance
    (
      GetObjectInstances().AddInstance( new GDT( 0.0))
    );

    //! @brief returns the default seed length
    //! @return the default seed length
    const size_t GDT::GetDefaultSeedLength()
    {
      return 3;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from a single distance cutoff and seed length
    //! @param DISTANCE_CUTOFF distance cutoff to be used
    //! @param SEED_LENGTH length of seed
    GDT::GDT
    (
      const double &DISTANCE_CUTOFF,
      const size_t SEED_LENGTH
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_SeedLength( GetDefaultSeedLength())
    {
    }

    //! @brief Clone function
    //! @return pointer to new GDT
    GDT *GDT::Clone() const
    {
      return new GDT( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GDT::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &GDT::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

    //! @brief get seed length
    //! @return seed length
    size_t GDT::GetSeedLength() const
    {
      return m_SeedLength;
    }

    //! @brief get distance cutoff
    //! @return distance cutoff
    double GDT::GetDistanceCutoff() const
    {
      return m_DistanceCutoff;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief collects the subset of coordinates specified by the given list of indices
    //! @param SUBSET indices of coordinates that define the subset of coordinates to be collected
    //! @param COORDINATES coordinates of interest
    //! @return vector of coordinates that correspond to the requested subset
    util::SiPtrVector< const linal::Vector3D> GDT::CollectCoordinatesSubset
    (
      const storage::List< size_t> &SUBSET,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    )
    {
      // initialize SiPtrVector to be returned
      util::SiPtrVector< const linal::Vector3D> coordinates_subset;
      coordinates_subset.AllocateMemory( SUBSET.GetSize());

      // iterate over the subset
      for
      (
        storage::List< size_t>::const_iterator index_itr( SUBSET.Begin()), index_itr_end( SUBSET.End());
        index_itr != index_itr_end; ++index_itr
      )
      {
        // insert the coordinates at the given index to the subset to be returned
        coordinates_subset.PushBack( COORDINATES( *index_itr));
      }

      // end
      return coordinates_subset;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates GDT between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return GDT between COORDINATES and REFERENCE_COORDINATES for a given distance cutoff
    double GDT::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // return the GDT value
      return CalculateGDTAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).First();
    }

    //! @brief calculates the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D GDT::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate GDT and return the superimposition
      return CalculateGDTAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).Second();
    }

    //! @brief calculates GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param DISTANCE_CUTOFF distance cutoff for the superimposition
    //! @return pair of GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    storage::Pair< double, math::TransformationMatrix3D> GDT::CalculateGDTAndSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // if the coordinates are empty are of different length
      if( COORDINATES.IsEmpty() || REFERENCE_COORDINATES.GetSize() != COORDINATES.GetSize())
      {
        return storage::Pair< double, math::TransformationMatrix3D>
        (
          util::GetUndefined< double>(), math::TransformationMatrix3D()
        );
      }

      // store number of points
      const size_t number_coordinates( COORDINATES.GetSize());

      // get the seed fragments
      storage::List< math::Range< size_t> > seed_fragments
      (
        GetSeedFragments( number_coordinates)
      );

      // get lcs as seeds also
      {
        const storage::List< math::Range< size_t> > lcs_seeds
        (
          LCS( m_DistanceCutoff).CalculateRanges( COORDINATES, REFERENCE_COORDINATES)
        );

        // if the lcs identified is larger than a normal seed add them
        if( !lcs_seeds.IsEmpty() && lcs_seeds.FirstElement().GetWidth() >= m_SeedLength)
        {
          seed_fragments.Append( lcs_seeds);
        }
      }

      // best indices and transformation
      storage::Pair< storage::List< size_t>, math::TransformationMatrix3D> best_indices_trans;

      // iterate over the seed fragments
      for
      (
        storage::List< math::Range< size_t> >::const_iterator itr( seed_fragments.Begin()), itr_end( seed_fragments.End());
           itr != itr_end
        && best_indices_trans.First().GetSize() != number_coordinates; // if the fragment includes all of the points meaning overall RMSD is less than the given cutoff
        ++itr
      )
      {
        // seed
        {
          const storage::Pair< storage::List< size_t>, math::TransformationMatrix3D> best_indices_trans_seed
          (
            IterativeSuperimpositionProcedure
            (
              COORDINATES, REFERENCE_COORDINATES,
              LCS::ConvertRangeToIndices( *itr)
            )
          );

          // if the longest fragment so far
          if( best_indices_trans_seed.First().GetSize() > best_indices_trans.First().GetSize())
          {
            // update it and best transformation
            best_indices_trans = best_indices_trans_seed;
          }
        }
      }

      // the longest fragment is as following
      BCL_MessageDbg
      (
        " the longest fragment and transformation is as following\n" + util::Format()( best_indices_trans)
      );

      // calculate the gdt value
      const double gdt( OptimalValue() * double( best_indices_trans.First().GetSize()) / double( number_coordinates));

      // return the pair of GDT value and the corresponding transformation matrix
      return storage::Pair< double, math::TransformationMatrix3D>( gdt, best_indices_trans.Second());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &GDT::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DistanceCutoff, ISTREAM);
      io::Serialize::Read( m_SeedLength    , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &GDT::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DistanceCutoff, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_SeedLength    , OSTREAM, 0);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create an average GDT object with given set of distance cutoffs
    //! @param DISTANCE_CUTOFFS set of distance cutoffs
    //! @param SEED_LENGTH length of seed
    //! @return Average object
    Average GDT::CreateAverageGDT
    (
      const storage::Set< double> &DISTANCE_CUTOFFS,
      const size_t SEED_LENGTH
    )
    {
      util::ShPtrVector< SuperimposeInterface> gdts;

      // for each distance cutoff given
      for( storage::Set< double>::const_iterator itr( DISTANCE_CUTOFFS.Begin()), itr_end( DISTANCE_CUTOFFS.End()); itr != itr_end; ++itr)
      {
        // construct a gdt measure
        gdts.PushBack( util::ShPtr< SuperimposeInterface>( new GDT( *itr, SEED_LENGTH)));
      }

      // end
      return Average( gdts);
    }

    //! @brief the iterative superimposition procedure (ISP)
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param START_INDICES the selected coordinates as seed, do not have to fullfill the CUTOFF
    //! @return the newly selected coordinates and its transformation
    storage::Pair< storage::List< size_t>, math::TransformationMatrix3D> GDT::IterativeSuperimpositionProcedure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
      const storage::List< size_t> &START_INDICES
    ) const
    {
      // a)
      const math::TransformationMatrix3D start_transform( CalculateTransformation( COORDINATES, REFERENCE_COORDINATES, START_INDICES));
      // b), c)
      storage::Pair< storage::List< size_t>, math::TransformationMatrix3D> best_indices_transform
      (
        CollectIndicesBelowCutoff( COORDINATES, REFERENCE_COORDINATES, start_transform),
        math::TransformationMatrix3D()
      );
      storage::List< size_t> &best_indices( best_indices_transform.First());

      // d)
      math::TransformationMatrix3D &best_transform( best_indices_transform.Second());
      best_transform = CalculateTransformation( COORDINATES, REFERENCE_COORDINATES, best_indices_transform.First());

      // nr_extensions
      size_t nr_extension( 0);

      while( true)
      {
        ++nr_extension;
        // b), c)
        storage::List< size_t> current_indices
        (
          CollectIndicesBelowCutoff( COORDINATES, REFERENCE_COORDINATES, best_transform)
        );
        // the identified indices are as good or worse than the best
        if( current_indices.GetSize() <= best_indices.GetSize())
        {
          break;
        }
        // d)
        best_indices.InternalData().swap( current_indices.InternalData());
        best_transform = CalculateTransformation( COORDINATES, REFERENCE_COORDINATES, best_indices);

        BCL_MessageDbg
        (
          "\textension # " + util::Format()( nr_extension) +
          " fragment of size " + util::Format()( best_indices.GetSize()) +
          " start indices:\n" + util::Format()( START_INDICES)
        );
      }

      // end
      return best_indices_transform;
    }

    //! @brief returns all possible seed ranges (not necessarily valid)
    //! @param NUMBER_COORDIANTES number of coordinates
    //! @return continuous seed fragments
    storage::List< math::Range< size_t> > GDT::GetSeedFragments
    (
      const size_t NUMBER_COORDINATES
    ) const
    {
      BCL_Assert( m_SeedLength >= 3, "The seed length should be at least 3! not: " + util::Format()( m_SeedLength));

      // initialize seed fragments list to return
      storage::List< math::Range< size_t> > seed_fragments;

      // check that there are enough coordinates
      if( NUMBER_COORDINATES < m_SeedLength)
      {
        return seed_fragments;
      }

      // iterate over beginning indices
      for( size_t index( 0), index_end( NUMBER_COORDINATES + 1 - m_SeedLength); index != index_end; ++index)
      {
        // create this seed
        const math::Range< size_t> this_seed( index, index + m_SeedLength - 1);

        // add it to seeds list
        seed_fragments.PushBack( this_seed);
      }

      // end
      return seed_fragments;
    }

    //! @brief check if two coordinates and the transformation will bring them close below cutoff
    //! @param COORDINATE coord to transform
    //! @param REFERENCE_COORDINATE coord to compare against
    //! @param TRANSFORMATION transformation to apply
    //! @param DISTANCE_CUTOFF_SQUARED the distance below which they are considered good
    //! @return true if the distance between the transformed coord and the reference is below cutoff
    bool GDT::IsBelowCutoff
    (
      const linal::Vector3D &COORDINATE,
      const linal::Vector3D &REFERENCE_COORDINATE,
      const math::TransformationMatrix3D &TRANSFORMATION,
      const double DISTANCE_CUTOFF_SQUARED
    )
    {
      // make a copy of the coordinate and transform
      linal::Vector3D coordinate_copy( COORDINATE);
      coordinate_copy.Transform( TRANSFORMATION);

      // calculate the distance squared between ref coordinate transformed coordinate
      // compare to distance cutoff
      return linal::SquareDistance( REFERENCE_COORDINATE, coordinate_copy) <= DISTANCE_CUTOFF_SQUARED;
    }

    //! @brief collect all indices below the cutoff
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param TRANSFORMATION transformation to apply before comparison
    //! @param list of coordinate indices that are below the cutoff for that transformation
    storage::List< size_t> GDT::CollectIndicesBelowCutoff
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATION
    ) const
    {
      // store number of points and number of seed fragments
      const size_t number_points( COORDINATES.GetSize());

      // store allowed distance cutoff and its squared value
      const double distance_cutoff_squared( m_DistanceCutoff * m_DistanceCutoff);

      // create a extended subset to return
      storage::List< size_t> indices;

      // now iterate over all residues
      for( size_t index( 0); index < number_points; ++index)
      {
        if( IsBelowCutoff( *COORDINATES( index), *REFERENCE_COORDINATES( index), TRANSFORMATION, distance_cutoff_squared))
        {
          // insert it into the new subset
          indices.PushBack( index);
        }
      } // end index

      // end
      return indices;
    }

    //! @brief transformation for a new set of indices
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @param INDICES the selected coordinates
    //! @return the transformation matrix to superimpose the coordinates selected by the indices
    math::TransformationMatrix3D GDT::CalculateTransformation
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES,
      const storage::List< size_t> &INDICES
    )
    {
      return RMSD::SuperimposeCoordinates
      (
        CollectCoordinatesSubset( INDICES, REFERENCE_COORDINATES),
        CollectCoordinatesSubset( INDICES, COORDINATES)
      );
    }

  } // namespace quality
} // namespace bcl
