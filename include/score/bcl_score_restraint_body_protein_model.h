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

#ifndef BCL_SCORE_RESTRAINT_BODY_PROTEIN_MODEL_H_
#define BCL_SCORE_RESTRAINT_BODY_PROTEIN_MODEL_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "coord/bcl_coord.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_body_assignment.h"
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintBodyProteinModel
    //! @brief class is for scoring the agreement of restraint::Body with the coord::Bodies
    //! of a protein model.
    //!
    //! @see @link example_score_restraint_body_protein_model.cpp @endlink
    //! @author alexanns
    //! @date 07/04/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintBodyProteinModel :
      public ProteinModel
    {
    /////////////
    // friends //
    /////////////

    private:

      //! the body restraints (i.e. the density rod restraints)
      util::ShPtr< util::ShPtrVector< restraint::Body> > m_Restraint;

      //! scores that will be used to score the agreement of the protein model with the body restraints "m_Restraint"
      BodyAssignment m_Score;

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintBodyProteinModel();

      //! @brief construct from a ShPtrList of Restraints and a score object
      //! @param RESTRAINT is the ShPtrList of Restraints which will be "m_Restraint"
      //! @param SCORE is the FunctionInterface object which will be used to score the restraints of "m_Restraint"
      RestraintBodyProteinModel
      (
        const util::ShPtr< util::ShPtrVector< restraint::Body> > &RESTRAINT,
        const BodyAssignment &SCORE
      );

      //! @brief Clone is the virtual copy constructor
      RestraintBodyProteinModel *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get scheme
      //! @return scheme
      const std::string &GetScheme() const;

      //!@brief GetRestraints returns a const reference to "m_Restraint"
      //!@return returns "m_Restraint"
      const util::ShPtr< util::ShPtrVector< restraint::Body> > &GetRestraints() const
      {
        return m_Restraint;
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief operator() which takes an ProteinModel for calculating its agreement with the restraint::Body
      //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the restraint::Body
      //! @return return a double which is the score of the agreement of the ProteinModel with the restraint::Body
      double operator()
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief ScoreCurrent gives the score the agreement of a ProteinModel with a Restraint
      //! @param RESTRAINT is the restraint which the proteinmodel is being score for agreement with
      //! @param PROTEIN_MODEL is the ProteinModel which is being scored for agreement with RESTRAINT
      //! @return returns a double which is the score of the agreement of PROTEIN_MODEL with RESTRAINT
      double ScoreCurrent
      (
        const restraint::Body &RESTRAINT,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class RestraintBodyProteinModel

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_BODY_PROTEIN_MODEL_H_
