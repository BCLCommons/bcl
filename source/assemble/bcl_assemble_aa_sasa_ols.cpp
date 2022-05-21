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
#include "assemble/bcl_assemble_aa_sasa_ols.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"
#include "coord/bcl_coord_sphere.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! @brief returns the default sphere radius
    //! @return the default sphere radius
    double AASasaOLS::GetDefaultSphereRadius()
    {
      // static integer to hold value
      static const double s_default_sphere_radius( 4.75);

      // end
      return s_default_sphere_radius;
    }

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &AASasaOLS::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "ols_sasa_all_chains_membrane.histograms");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &AASasaOLS::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "sasaols");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param SCHEME scheme to be used
    AASasaOLS::AASasaOLS( const std::string &SCHEME) :
      m_Scheme( SCHEME),
      m_SphereRadiusThreshold( 0, 2 * GetDefaultSphereRadius()),
      m_MininalSequenceSeparation( s_DefaultMinimalSequenceSeparation)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AASasaOLS object copied from this one
    AASasaOLS *AASasaOLS::Clone() const
    {
      return new AASasaOLS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASasaOLS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AASasaOLS::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief min and max exposure measure
    //! @return range in which exposure can be
    const math::Range< double> &AASasaOLS::GetRange() const
    {
      static const math::Range< double> s_range( 0, 1);
      return s_range;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate SASA for a given amino acid and its AANeighborList
    //! @return SASA for a given amino acid and its AANeighborList
    double AASasaOLS::operator()( const AANeighborList &AA_NEIGHBOR_LIST) const
    {
      // radius
      const double radius( m_SphereRadiusThreshold.GetMax() / 2);

      // check that the given neighbor list has the proper parameters
      BCL_Assert
      (
        AA_NEIGHBOR_LIST.GetDistanceCutoff() == 2 * radius,
        "given neighbor list does not have the proper distance cutoff: " +
        util::Format()( AA_NEIGHBOR_LIST.GetDistanceCutoff()) + " == " + util::Format()( 2 * radius)
      );

      if( !AA_NEIGHBOR_LIST.GetCenterAminoAcid()->GetFirstSidechainAtom().GetCoordinates().IsDefined())
      {
        return util::GetUndefined< double>();
      }

      // sphere for center aa
      const coord::Sphere center_sphere
      (
        AA_NEIGHBOR_LIST.GetCenterAminoAcid()->GetFirstSidechainAtom().GetCoordinates(), radius
      );

      // create a sphere for all neighbor amino acids
      storage::List< coord::Sphere> spheres;

      // iterate over all neighbors
      for
      (
        AANeighborList::const_iterator itr( AA_NEIGHBOR_LIST.Begin()), itr_end( AA_NEIGHBOR_LIST.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->First()->GetFirstSidechainAtom().GetCoordinates().IsDefined())
        {
          spheres.PushBack
          (
            coord::Sphere( itr->First()->GetFirstSidechainAtom().GetCoordinates(), radius)
          );
        }
      }

      // calculate the free surface area fraction
      return center_sphere.FreeSurfaceAreaFraction( util::ConvertToConstSiPtrList< const coord::Sphere>( spheres.Begin(), spheres.End()), 100);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AASasaOLS::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme               , ISTREAM);
      io::Serialize::Read( m_SphereRadiusThreshold, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @return returns the output stream
    std::ostream &AASasaOLS::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme               , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SphereRadiusThreshold, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
