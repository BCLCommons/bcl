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
#include "score/bcl_score_protein_model_completeness.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelCompleteness::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelCompleteness())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from default values
    //! @param IGNORE_TERM_LOOPS ignore terminal loops for the calculation of the completeness
    //! @param SCHEME scheme of the score
    ProteinModelCompleteness::ProteinModelCompleteness( bool IGNORE_TERM_LOOPS, const std::string &SCHEME) :
      m_IgnoreTermLoops( IGNORE_TERM_LOOPS),
      m_Scheme( SCHEME)
    {
    }

    //! @brief returns a pointer to a new ProteinModelCompleteness
    //! @return pointer to a new ProteinModelCompleteness
    ProteinModelCompleteness *ProteinModelCompleteness::Clone() const
    {
      return new ProteinModelCompleteness( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &ProteinModelCompleteness::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the scheme of this score
    //! @return the scheme of this score
    const std::string &ProteinModelCompleteness::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief returns the default scheme of this score
    //! @return the default scheme of this score
    const std::string &ProteinModelCompleteness::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "coord_complete");
      return s_default_scheme;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProteinModelCompleteness::GetAlias() const
    {
      static const std::string s_name( "ProteinModelCompleteness");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelCompleteness::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Scores the completeness of a protein model.");
      serializer.AddInitializer
      (
        "ignore term loops",
        "ignore terminal loops for the completeness computation",
        io::Serialization::GetAgent( &m_IgnoreTermLoops),
        "False"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief scores the completeness of the given protein model
    //! @detail the completeness is calculated by computing the ratio of residues with defined coordinates and the
    //! total number of residues in the sequence.
    //! @param PROTEIN_MODEL protein model for which to compute the completeness
    //! @return completeness of the given protein model with -1.0 being complete and 0.0 being empty
    double ProteinModelCompleteness::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // compute the total number of residues and the residues with defined coordinates
      const size_t num_residues( GetNumResidues( PROTEIN_MODEL));
      const size_t num_residues_defined( GetNumResiduesDefined( PROTEIN_MODEL));

      // compute the fraction only if there are residues in the model
      const double score( num_residues != 0 ? -( double) num_residues_defined / num_residues : 0.0);

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief computes the number of residues in the given chain
    //! @param CHAIN chain for which to compute the number of residues
    //! @return number of residues with in the given chain
    size_t ProteinModelCompleteness::GetNumResidues( const assemble::Chain &CHAIN) const
    {
      // get the total number of residues in the chain
      size_t num_residues( CHAIN.GetNumberAAs());

      // subtract the residues in terminal loops if required
      if( m_IgnoreTermLoops)
      {
        // skip if the first SSE is not a loop, otherwise subtract the number of residues in the loop
        const assemble::SSE &first_sse( **CHAIN.GetData().Begin());
        if( first_sse.GetType() == biol::GetSSTypes().COIL)
        {
          num_residues -= first_sse.GetSize();
        }

        // skip if the last SSE is not a loop, otherwise subtract the number of residues in the loop
        const assemble::SSE &last_sse( **--CHAIN.GetData().End());
        if( last_sse.GetType() == biol::GetSSTypes().COIL)
        {
          num_residues -= last_sse.GetSize();
        }
      }

      return num_residues;
    }

    //! @brief computes the number of residues in the given protein model
    //! @param MODEL protein model for which to compute the number of residues
    //! @return number of residues in the given protein model
    size_t ProteinModelCompleteness::GetNumResidues( const assemble::ProteinModel &MODEL) const
    {
      // sum up the number of residues in each chain
      size_t num_residues( 0);
      for
      (
        auto chain_it( MODEL.GetChains().Begin()), chain_it_end( MODEL.GetChains().End());
        chain_it != chain_it_end;
        ++chain_it
      )
      {
        num_residues += GetNumResidues( **chain_it);
      }

      return num_residues;
    }

    //! @brief computes the number of residues with defined backbone coordinates in a given SSE
    //! @param SSE SSE for which to compute the number of residues with defined backbone coordinates
    //! @return number of residues with defined backbone coordinates in the given SSE
    size_t ProteinModelCompleteness::GetNumResiduesDefined( const assemble::SSE &SSE) const
    {
      // get the residues in the SSE
      const util::SiPtrVector< const biol::AABase> &residues( SSE.GetData());

      // count the residues with defined backbone coordinates
      size_t num_residues_defined( 0);
      for( auto res_it( residues.Begin()), res_it_end( residues.End()); res_it != res_it_end; ++res_it)
      {
        if( ( **res_it).HasDefinedCoordinates())
        {
          ++num_residues_defined;
        }
      }

      return num_residues_defined;
    }

    //! @brief computes the number of residues with defined backbone coordinates in a given chain
    //! @param CHAIN chain for which to compute the number of residues with defined backbone coordinates
    //! @return number of residues with defined backbone coordinates in the given chain
    size_t ProteinModelCompleteness::GetNumResiduesDefined( const assemble::Chain &CHAIN) const
    {
      // get the residues in the chain
      const util::SiPtrVector< const biol::AABase> &residues( CHAIN.GetAminoAcids());

      // count the residues with defined backbone coordinates
      size_t num_residues_defined( 0);
      for( auto res_it( residues.Begin()), res_it_end( residues.End()); res_it != res_it_end; ++res_it)
      {
        if( ( **res_it).HasDefinedCoordinates())
        {
          ++num_residues_defined;
        }
      }

      // subtract the defined residues in terminal loops if required
      if( m_IgnoreTermLoops)
      {
        // skip if the first SSE is not a loop, otherwise subtract the number of residues in the loop
        const assemble::SSE &first_sse( *CHAIN.GetSSEs().FirstElement());
        if( first_sse.GetType() == biol::GetSSTypes().COIL)
        {
          num_residues_defined -= GetNumResiduesDefined( first_sse);
        }

        // skip if the last SSE is not a loop, otherwise subtract the number of residues in the loop
        const assemble::SSE &last_sse( *CHAIN.GetSSEs().LastElement());
        if( last_sse.GetType() == biol::GetSSTypes().COIL)
        {
          num_residues_defined -= GetNumResiduesDefined( last_sse);
        }
      }

      return num_residues_defined;
    }

    //! @brief computes the number of residues with defined backbone coordinates in a given protein model
    //! @param MODEL protein model for which to compute the number of residues with defined backbone coordinates
    //! @return number of residues with defined backbone coordinates in the given protein model
    size_t ProteinModelCompleteness::GetNumResiduesDefined( const assemble::ProteinModel &MODEL) const
    {
      // sum up the number of defined residues in each chain
      size_t num_residues_defined( 0);
      const util::ShPtrVector< assemble::Chain> chains( MODEL.GetChains());
      for( auto chain_it( chains.Begin()), chain_it_end( chains.End()); chain_it != chain_it_end; ++chain_it)
      {
        num_residues_defined += GetNumResiduesDefined( **chain_it);
      }

      return num_residues_defined;
    }

  } // namespace score
} // namespace bcl
