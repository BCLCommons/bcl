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
#include "score/bcl_score_protein_model_defined_loops.h"

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
    const util::SiPtr< const util::ObjectInterface> ProteinModelDefinedLoops::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelDefinedLoops())
    );
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from default values
    //! @param IGNORE_TERM_LOOPS ignore terminal loops for the calculation of the completeness
    //! @param SCHEME scheme of the score
    ProteinModelDefinedLoops::ProteinModelDefinedLoops( bool IGNORE_TERM_LOOPS, const std::string &SCHEME) :
      m_IgnoreTermLoops( IGNORE_TERM_LOOPS),
      m_Scheme( SCHEME)
    {
    }

    //! @brief returns a pointer to a new ProteinModelDefinedLoops
    //! @return pointer to a new ProteinModelDefinedLoops
    ProteinModelDefinedLoops *ProteinModelDefinedLoops::Clone() const
    {
      return new ProteinModelDefinedLoops( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &ProteinModelDefinedLoops::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the scheme of this score
    //! @return the scheme of this score
    const std::string &ProteinModelDefinedLoops::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief returns the default scheme of this score
    //! @return the default scheme of this score
    const std::string &ProteinModelDefinedLoops::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "loops_defined");
      return s_default_scheme;
    }
    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &ProteinModelDefinedLoops::GetAlias() const
    {
      static const std::string s_name( "ProteinModelDefinedLoops");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelDefinedLoops::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores the loop completeness of a protein model");
      serializer.AddInitializer
      (
        "ignore term loops",
        "ignore terminal loops for the calculation of the completeness",
        io::Serialization::GetAgent( &m_IgnoreTermLoops),
        "false"
      );
      return serializer;
    }
  ////////////////
  // operations //
  ////////////////

    //! @brief scores the loop completeness of the given protein model
    //! @detail the loop completeness is calculated by computing the number of loops with fully defined
    //! backbone coordinates.
    //! @param PROTEIN_MODEL protein model for which to compute the loop completeness
    //! @return loop completeness of the given protein model with -1 * number of defined loops
    double ProteinModelDefinedLoops::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // get the loops in the protein model
      const util::SiPtrVector< const assemble::SSE> loops( PROTEIN_MODEL.GetSSEs( biol::GetSSTypes().COIL));

      // determine the number of loops with undefined coordinates
      size_t undefined_loops( 0);
      for( auto loop_it( loops.Begin()), loop_it_end( loops.End()); loop_it != loop_it_end; ++loop_it)
      {
        // loop to be scored
        const assemble::SSE &loop( **loop_it);

        // check if current loop is terminal and ignore it if necessary
        const biol::AABase &first_res( *loop.GetFirstMember());
        const biol::AABase &last_res( *loop.GetLastMember());
        const size_t seq_length
        (
          PROTEIN_MODEL.GetChain( loop.GetChainID())->GetSequence()->GetMembers().GetSize()
        );
        if( !m_IgnoreTermLoops || ( first_res.GetSeqID() != 1 && ( size_t) last_res.GetSeqID() != seq_length))
        {
          if( !IsDefined( loop))
          {
            ++undefined_loops;
          }
        }
      }

      return undefined_loops;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from an input stream
    //! @param ISTREAM input stream to read members from
    //! @return the input stream
    std::istream &ProteinModelDefinedLoops::Read( std::istream &ISTREAM)
    {
      // read members from input stream
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_IgnoreTermLoops, ISTREAM);

      return ISTREAM;
    }

    //! @brief writes members into an output stream
    //! @param OSTREAM output stream to write members into
    //! @INDENT number of indentations to use
    //! @return the output stream
    std::ostream &ProteinModelDefinedLoops::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members into output stream
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_IgnoreTermLoops, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief determines if all coordinates in the given SSE are defined
    //! @param SSE SSE for which to determine if all coordinates are defined
    //! @return true, if all coordinates of the given SSE are defined
    bool ProteinModelDefinedLoops::IsDefined( const assemble::SSE &SSE) const
    {
      // get the residues in the SSE
      const util::SiPtrVector< const biol::AABase> &residues( SSE.GetData());

      // count the residues with defined backbone coordinates
      for( auto res_it( residues.Begin()), res_it_end( residues.End()); res_it != res_it_end; ++res_it)
      {
        if( !( **res_it).HasDefinedCoordinates())
        {
          return false;
        }
      }

      return true;
    }

  } // namespace score
} // namespace bcl
