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
#include "score/bcl_score_restraint_noe_attraction.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // map of histogram data
    storage::Map
    <
      storage::Pair< biol::AtomType, size_t>,
      RestraintAtomAttraction
    > RestraintNoeAttraction::s_HistogramData;

    // initialize default score
    const double RestraintNoeAttraction::s_DefaultScore( 0.0);

    const util::SiPtr< const util::ObjectInterface> RestraintNoeAttraction::s_Instance
    (
      util::Enumerated< RestraintAtomDistanceAssignment>::AddInstance( new RestraintNoeAttraction())
    );

  //////////
  // data //
  //////////

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &RestraintNoeAttraction::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "noe_knowledge_based.histogram");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &RestraintNoeAttraction::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "atom_attraction_noe");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual copy constructor
    RestraintNoeAttraction *RestraintNoeAttraction::Clone() const
    {
      return new RestraintNoeAttraction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RestraintNoeAttraction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the histograms
    //! @return the histograms
    const storage::Map
    <
      storage::Pair< biol::AtomType, size_t>,
      RestraintAtomAttraction
    > &RestraintNoeAttraction::GetNOEHistogram()
    {
      // read in the data if it is empty
      if( s_HistogramData.IsEmpty())
      {
        s_HistogramData = ReadNOEHistogram();
      }

      // end
      return s_HistogramData;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator scores protein model
    //! @param RESTRAINT restraint to be scored
    //! @return score
    double RestraintNoeAttraction::operator()( const restraint::AtomDistanceAssignment &RESTRAINT) const
    {
      // calculate the distance
      const double cb_distance( RESTRAINT.CalculateAtomDistance());

      // if the calculated distance is undefined
      if( !util::IsDefined( cb_distance))
      {
        // return default score
        return s_DefaultScore;
      }

      // get the bond distance
      const size_t bond_distance( GetTotalBondsFromCB( RESTRAINT));

      // if the bond distance is not defined
      if( !util::IsDefined( bond_distance))
      {
        // return default score
        return s_DefaultScore;
      }

      // return the value calculated from the appropriate spline
      return GetNOEHistogram().GetValue
      (
        storage::Pair< biol::AtomType, size_t>( GetAtomType( RESTRAINT), bond_distance)
      )( RESTRAINT);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief reads in the histograms with NOE distance statistics from the file
    //! @return the histograms
    storage::Map
    <
      storage::Pair< biol::AtomType, size_t>,
      RestraintAtomAttraction
    > RestraintNoeAttraction::ReadNOEHistogram()
    {
      // create IFStream "read"
      io::IFStream read;

      // open "read" and bind it to the histogram file containing SL-CB distances
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( GetDefaultHistogramFilename()));

      // initialize atom type
      biol::AtomType atom_type;

      // initialize map
      storage::Map
      <
        storage::Pair< biol::AtomType, size_t>,
        RestraintAtomAttraction
      > histograms;

      // iterate through the entire histogram
      while( read >> atom_type && !read.eof())
      {
        // break if undefined atom_type (end of file)
        if( atom_type == biol::GetAtomTypes().e_Undefined)
        {
          break;
        }

        // create an empty size_t for storing the bond distances
        size_t bond_distance;

        // read histogram
        math::Histogram histogram;
        read >> bond_distance >> histogram;

        // calculate spline and store in map
        histograms[ storage::Pair< biol::AtomType, size_t>( atom_type, bond_distance)] =
          RestraintAtomAttraction
          (
            RestraintAtomAttraction::GetDefaultDepthRange(),
            RestraintAtomAttraction::GetDefaultLeftEndWell( histogram),
            RestraintAtomAttraction::GetDefaultTransitionWidth(),
            true
          );
      }

      // close and clear read stream
      io::File::CloseClearFStream( read);

      // end
      return histograms;
    }

  } // namespace score
} // namespace bcl
