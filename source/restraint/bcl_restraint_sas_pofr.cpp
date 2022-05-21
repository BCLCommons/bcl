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
#include "restraint/bcl_restraint_sas_pofr.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_ca_cb.h"
#include "biol/bcl_biol_atom_group_types.h"
#include "biol/bcl_biol_sasa_data.h"
#include "fold/bcl_fold_add_parabolic_loops.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_sum_function.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_density.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SasPofR::s_Instance
    (
      util::Enumerated< SasPofRInterface>::AddInstance( new SasPofR())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new SasPofR function
    SasPofR *SasPofR::Clone() const
    {
      return new SasPofR( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SasPofR::GetAlias() const
    {
      static const std::string s_Name( "SasPofR");
      return s_Name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasPofR::GetClassIdentifier() const
    {
      // Get BCL Standardized Class name
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Get the positions for all the atoms in the protein model
    //! @param MODEL the protein model of interest
    //! @return a vector containing positions
    storage::Vector< linal::Vector3D> SasPofR::GetAtoms( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      BCL_MessageDbg( " Inside Get all atom coordinates from the protein model: ");

      // initialize vector to hold atom coordinates
      storage::Vector< linal::Vector3D> all_atoms;

      size_t atom_number( 0);

      // iterate over the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
         chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // get all the structured sses
        const util::SiPtrVector< const assemble::SSE> structured_sses
        (
          ( *chain_itr)->GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
        );

        // if structured_sses are empty
        if( structured_sses.IsEmpty())
        {
          // warn user and return empty vector
          BCL_Message( util::Message::e_Standard, "No structured SSEs found in protein model");
          continue;
        }

        // iterate over all the SSEs including coil (depending on min_sse_size parameter)
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
           sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // iterate over the Amino Acids in a given SSE
          for
          (
            util::ShPtrVector< biol::AABase>::const_iterator aa_itr( ( *sse_itr)->GetData().Begin()),
             aa_itr_end( ( *sse_itr)->GetData().End());
            aa_itr != aa_itr_end; ++aa_itr
          )
          {
            std::string amino_acid( ( *aa_itr)->GetType()->GetName());
            BCL_MessageDbg( "Residue: " + util::Format()( amino_acid));

            // get the atom types for this residue
            // Be sure to get the Heavy Atom Types ( all atoms except hydrogen).  This set contains all of the residues
            // That a given amino acid contains.  Not the atoms actually present in the input PDB file
            const storage::Set< biol::AtomType> &complete_atom_types( ( *aa_itr)->GetType()->GetAllowedHeavyAtomTypes());

            BCL_MessageDbg( "AtomTypes: " + util::Format()( complete_atom_types));

            // iterate over the atoms in a given residue
            for
            (
              storage::Set< biol::AtomType>::const_iterator atom_itr( complete_atom_types.Begin()),
               atom_itr_end( complete_atom_types.End());
              atom_itr != atom_itr_end; ++atom_itr
            )
            {
              atom_number++;

              // get the atom coordinates
              const linal::Vector3D &atom_coords( ( *aa_itr)->GetAtom( *atom_itr).GetCoordinates());

              BCL_MessageDbg( "Atom Coords: " + util::Format()( atom_coords));

              if( atom_coords.IsDefined())
              {
                std::string atom( ( *atom_itr)->GetName());
              }

              // if the coordinates are defined
              if( atom_coords.IsDefined())
              {
                // pushback into vector
                all_atoms.PushBack( atom_coords);
              }
            } // Atom Type Iteration
          } // AA iteration
        } // SSE iteration
      } // chain iteration
      return all_atoms;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a proteinModel as input and calculates an intensity using the debye formula
    //! @param PROTEIN_MODEL
    //! @return the intensity for a given q value
    SasExperimentalAndCalculatedDensity SasPofR::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // This holds the Euclidean coordinates of the atoms in the model
      storage::Vector< linal::Vector3D> all_atoms( GetAtoms( PROTEIN_MODEL));

      BCL_MessageDbg( " All atoms: " + util::Format()( all_atoms));

      double bin_size( this->GetExperimentalDensity()->GetBinSize());
      size_t bin_number( this->GetExperimentalDensity()->GetBinNumber());

      // construct Histogram from starting bin, bin size and number of bins
      math::Histogram distance_histogram( 0.0, bin_size, bin_number);

      size_t row( 0);
      double dmax( 0);

      // Compute the upper triangle of pairwise distances and increment the historgram
      for
      (
        storage::Vector< linal::Vector3D>::const_iterator atom_itr_a( all_atoms.Begin()), atom_itr_end( all_atoms.End());
        atom_itr_a != atom_itr_end; ++atom_itr_a, ++row
      )
      {
        size_t col( row + 1);
        for
        (
          storage::Vector< linal::Vector3D>::const_iterator atom_itr_b( atom_itr_a + 1);
          atom_itr_b != atom_itr_end; ++atom_itr_b, ++col
        )
        {

          double distance = linal::Distance( *atom_itr_a, *atom_itr_b);
          distance_histogram.PushBack( distance);
          if( dmax < distance)
          {
            dmax = distance;
          }

          // Take care of the lower half of the matrix
          distance_histogram.PushBack( distance);
        }
      }

      // create object to hold calculated data
      // The computed dmax is set here, the experimental dmax is read in at object creation
      SasDensityData calculated_data( distance_histogram, dmax);

      SasExperimentalAndCalculatedDensity pofr_data( *this->GetExperimentalDensity(), calculated_data);

      BCL_MessageDbg( "data: " + util::Format()( pofr_data));

      double cal_dmax( pofr_data.GetCalculatedDensity().GetDmax());
      double exp_dmax( pofr_data.GetExperimentalDensity().GetDmax());

      BCL_MessageDbg( "cal_dmax: " + util::Format()( cal_dmax));
      BCL_MessageDbg( "exp_dmax: " + util::Format()( exp_dmax));

      double cal_hmax( pofr_data.GetCalculatedDensity().GetHmax());
      double exp_hmax( pofr_data.GetExperimentalDensity().GetHmax());

      BCL_MessageDbg( "cal_hmax: " + util::Format()( cal_hmax));
      BCL_MessageDbg( "exp_hmax: " + util::Format()( exp_hmax));

      return pofr_data;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SasPofR::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "performs pairwise atomic distance calculation"
      );

      return parameters;
    }

  } // namespace restraint
} // namespace bcl
