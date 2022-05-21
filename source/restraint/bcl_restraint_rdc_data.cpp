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
#include "restraint/bcl_restraint_rdc_data.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_score_weight_set.h"
#include "nmr/bcl_nmr_residual_dipolar_coupling_least_square_deviation.h"
#include "nmr/bcl_nmr_star_rdc_handler.h"
#include "score/bcl_score_residual_dipolar_coupling_q_value.h"
#include "score/bcl_score_restraint_residual_dipolar_coupling.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // initialize score
    fold::Score RDCData::e_ScoreRDCRestraint( fold::GetScores().e_Undefined);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RDCData::s_Instance
    (
      util::Enumerated< Interface>::AddInstance( new RDCData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RDCData::RDCData() :
      m_Restraints(),
      m_Handler( GetDefaultHandler())
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param HANDLER the handler that will be used to read and write the restraints
    RDCData::RDCData( const HandlerBase< RDC> &HANDLER) :
      m_Restraints(),
      m_Handler( HANDLER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RDCData
    RDCData *RDCData::Clone() const
    {
      return new RDCData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &RDCData::GetAlias() const
    {
      static const std::string s_name( "RDC");
      return s_name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RDCData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the default file extension
    //! @return the default file extension
    const std::string &RDCData::GetDefaultExtension() const
    {
      static const std::string s_extension( ".rdc_star");
      return s_extension;
    }

    //! @brief get the default restraint handler
    //! @return the default restraint handler
    const nmr::StarRDCHandler &RDCData::GetDefaultHandler()
    {
      static const nmr::StarRDCHandler s_handler( ".rdc_star");
      return s_handler;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void RDCData::InitializeScores()
    {
      if( !e_ScoreRDCRestraint.IsDefined())
      {
        if( !m_Restraints.IsDefined() && m_Handler.IsDefined() && m_Handler->Exists())
        {
          m_Restraints = util::CloneToShPtr( m_Handler->ReadRestraintsFromFile());
        }
        else
        {
          *m_Restraints = m_Handler->ReadRestraintsFromFile();
        }

        e_ScoreRDCRestraint = fold::GetScores().AddScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::RestraintResidualDipolarCoupling
            (
              m_Restraints,
              nmr::ResidualDipolarCouplingLeastSquareDeviation(),
              score::ResidualDipolarCouplingQValue()
            )
          )
        );
      }
    }

    //! @brief sets the weights of scores in a weight set
    //! @param SCORE_WEIGHT_SET the score weight set that will be modified
    void RDCData::ModifyScoreWeightSet( fold::ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.SetWeight( e_ScoreRDCRestraint, 5);
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void RDCData::MergeAndModifyMutateTree( fold::MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< fold::MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }
 
  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RDCData::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription
      (
        "Residual dipolar coupling restraints"
      );
      serial.AddInitializer
      (
        "",
        "Handler for reading RDCs",
        io::Serialization::GetAgent( &m_Handler),
        util::Implementation< HandlerBase< RDC> >( GetDefaultHandler()).GetString()
      );
      return serial;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RDCData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Restraints, ISTREAM);
      io::Serialize::Read( m_Handler, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RDCData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Handler, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
