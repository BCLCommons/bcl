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
#include "density/bcl_density_protein_agreement_likelihood.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_mask_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////
  // data //
  //////////

    //! @brief low res masking distance for Mask3D using CA only
    const double ProteinAgreementLikelihood::s_LowResolutionMaskingDistance( 8.0);

    //! @brief high res masking distance for Mask3D using all side chain atoms
    const double ProteinAgreementLikelihood::s_HighResolutionMaskingDistance( 5.0);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinAgreementLikelihood::s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinAgreementLikelihood::ProteinAgreementLikelihood() :
      m_LogLikelihood(),
      m_HighResolution( false),
      m_AtomTypes( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA)),
      m_DensityMap(),
      m_Simulator()
    {
      SetScheme();
    }

    //! @brief constructor from atom types
    //! @param HIGH_RESOLUTION use all atoms of residue plus neighboring residues, otherwise just CA
    //! @param ATOM_TYPES atom types to consider
    //! @param MEAN_CCC_RESOLUTION_FUNCTION
    //! @param SD_CCC_RESOLUTION_FUNCTION
    ProteinAgreementLikelihood::ProteinAgreementLikelihood
    (
      const bool HIGH_RESOLUTION,
      const storage::Set< biol::AtomType> &ATOM_TYPES,
      const math::FunctionInterfaceSerializable< double, double> &MEAN_CCC_RESOLUTION_FUNCTION,
      const math::FunctionInterfaceSerializable< double, double> &SD_CCC_RESOLUTION_FUNCTION
    ) :
      m_HighResolution( HIGH_RESOLUTION),
      m_AtomTypes( ATOM_TYPES),
      m_MeanFit( MEAN_CCC_RESOLUTION_FUNCTION.Clone()),
      m_SdFit( SD_CCC_RESOLUTION_FUNCTION.Clone())
    {
      SetScheme();
    }

    //! @brief Clone function
    //! @return pointer to new ProteinAgreementLikelihood
    ProteinAgreementLikelihood *ProteinAgreementLikelihood::Clone() const
    {
      return new ProteinAgreementLikelihood( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinAgreementLikelihood::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinAgreementLikelihood::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief access to the density simulator
    //! @return the density simulator used
    const util::ShPtr< SimulateInterface> &ProteinAgreementLikelihood::GetSimulator() const
    {
      return m_Simulator;
    }

    //! @brief set the simulator used for density agreement, e.g. for simulating density from the protein model
    //! @param SP_SIMULATOR ShPtr to SimulatInterface
    void ProteinAgreementLikelihood::SetSimulator( const util::ShPtr< SimulateInterface> &SP_SIMULATOR)
    {
      m_Simulator = SP_SIMULATOR;
      const double resolution( SP_SIMULATOR->GetResolution());
      m_LogLikelihood = math::LogLikelihood( m_MeanFit->operator ()( resolution), m_SdFit->operator ()( resolution));
    }

    //! @brief access to the density used for agreement calculation
    //! @return SiPtr to the density
    const util::SiPtr< const Map> &ProteinAgreementLikelihood::GetDensity() const
    {
      return m_DensityMap;
    }

    //! @brief set the density used for agreement calculation
    //! @param SP_DENSITY SiPtr to the density map
    void ProteinAgreementLikelihood::SetDensityMap( const util::SiPtr< const Map> &SP_DENSITY)
    {
      m_DensityMap = SP_DENSITY;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the standard Deviation between the protein model and the given density
    //! @param PROTEIN_MODEL protein of interest
    //! @return score for the density likelihood
    double ProteinAgreementLikelihood::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // high resolution score
      if( m_HighResolution)
      {
        return HighResolutionScore( PROTEIN_MODEL);
      }
      // low resolution score
      else
      {
        return LowResolutionScore( PROTEIN_MODEL);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinAgreementLikelihood::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DensityMap    , ISTREAM);
      io::Serialize::Read( m_HighResolution, ISTREAM);
      io::Serialize::Read( m_AtomTypes     , ISTREAM);
      io::Serialize::Read( m_LogLikelihood , ISTREAM);
      io::Serialize::Read( m_Simulator     , ISTREAM);

      // set scheme according to atom types
      SetScheme();

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinAgreementLikelihood::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_DensityMap    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HighResolution, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomTypes     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LogLikelihood , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Simulator     , OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate low resolution score, only considering CA atoms
    //! @param PROTEIN_MODEL protein of interest
    //! @return score for the density likelihood
    double ProteinAgreementLikelihood::LowResolutionScore( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // simulate density for given PROTEIN_MODEL
      const Map density_sim( m_Simulator->operator()( PROTEIN_MODEL.GetAtoms( m_AtomTypes)));

      // get all residues
      const util::SiPtrVector< const biol::AABase> residues( PROTEIN_MODEL.GetAminoAcids());

      // score
      double score( 0.0);

      // iterate over all residues
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator aa_itr( residues.Begin()), aa_itr_end( residues.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // current res coordinates
        const util::SiPtrVector< const linal::Vector3D> coords( ( *aa_itr)->GetAtomCoordinates( m_AtomTypes));

        // skip undefined, as for GLY
        if( !coord::AreDefinedCoordinates( coords))
        {
          continue;
        }

        // generate mask with coords
        const Mask3d current_mask
        (
          coords, s_LowResolutionMaskingDistance, density_sim.GetCellWidth(), density_sim.GetOrigin()
        );

        // calculate cross correlation over mask
        const double mask_ccc( current_mask.CrossCorrelationCoefficient( *m_DensityMap, density_sim));
        if( util::IsDefined( mask_ccc))
        {
          // evaluate log likelihood function and add to score
          score += m_LogLikelihood( mask_ccc);
        }
      }

      // end
      return score;
    }

    //! @brief calculate high resolution score, considering all atoms of each amino acid sidechain plus the two
    //! neighboring residues
    //! @param PROTEIN_MODEL protein of interest
    //! @return score for the density likelihood
    double ProteinAgreementLikelihood::HighResolutionScore( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // generate new protein model with side chains from PROTEIN_MODEL
      const assemble::ProteinModel model
      (
        *biol::AASideChainFactory( false, true).ProteinModelWithSideChains( PROTEIN_MODEL)
      );

      // get all residues
      const util::SiPtrVector< const biol::AABase> residues( model.GetAminoAcids());

      // score
      double score( 0.0);

      // need at least 3 amino acids for high res score
      if( residues.GetSize() < 3)
      {
        return score;
      }

      // simulate density for given PROTEIN_MODEL
      const Map density_sim( m_Simulator->operator()( model.GetAtoms()));

      // iterate over all residues
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          aa_itr_left( residues.Begin()), aa_itr_current( aa_itr_left + 1), aa_itr_right( aa_itr_current + 1),
          aa_itr_end( residues.End());
        aa_itr_right != aa_itr_end;
        ++aa_itr_left, ++aa_itr_current, ++aa_itr_right
      )
      {
        const biol::AABase &aa_left( **aa_itr_left);
        const biol::AABase &aa_current( **aa_itr_current);
        const biol::AABase &aa_right( **aa_itr_right);

        // sequence separation
        const size_t seq_separation_lc( biol::SequenceSeparation( aa_left, aa_current));
        const size_t seq_separation_cr( biol::SequenceSeparation( aa_current, aa_right));

        // check that all amino acids are from the same chain and are neighboring amino acids
        if( seq_separation_lc != 0 || seq_separation_cr != 0)
        {
          continue;
        }

        util::SiPtrVector< const linal::Vector3D> coords;
        coords.Append( aa_left.GetAtomCoordinates());
        coords.Append( aa_current.GetAtomCoordinates());
        coords.Append( aa_right.GetAtomCoordinates());
        const Mask3d current_mask
        (
          coords, s_HighResolutionMaskingDistance, density_sim.GetCellWidth(), density_sim.GetOrigin()
        );

        // calculate cross correlation over mask
        const double mask_ccc( current_mask.CrossCorrelationCoefficient( *m_DensityMap, density_sim));
        if( util::IsDefined( mask_ccc))
        {
          // evaluate log likelihood function and add to score
          score += m_LogLikelihood( mask_ccc);
        }
      }

      // end
      return score;
    }

    //! @brief string for the default scheme
    const std::string &ProteinAgreementLikelihood::GetDefaultScheme() const
    {
      // initialize static scheme
      static const std::string s_scheme( "Likelihood");

      // end
      return s_scheme;
    }

    //! @brief set scheme
    //! @details set scheme from the default scheme and atom types
    void ProteinAgreementLikelihood::SetScheme()
    {
      m_Scheme = GetDefaultScheme();
      for
      (
        storage::Set< biol::AtomType>::const_iterator itr( m_AtomTypes.Begin()), itr_end( m_AtomTypes.End());
        itr != itr_end;
        ++itr
      )
      {
        m_Scheme += ( *itr)->GetName();
      }
    }

  } // namespace density
} // namespace bcl
