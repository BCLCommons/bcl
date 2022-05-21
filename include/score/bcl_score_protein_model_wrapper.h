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

#ifndef BCL_SCORE_PROTEIN_MODEL_WRAPPER_H_
#define BCL_SCORE_PROTEIN_MODEL_WRAPPER_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelWrapper
    //! @brief Protein model score from math::FunctionInterfaceSerializable< assemble::ProteinModel, double>
    //!
    //! @remarks example unnecessary
    //! @author weinerbe
    //! @date Jul 7, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelWrapper :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! function to wrap
      util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, double> > m_Function;

      //! score type
      Type m_ScoreType;

      //! readable scheme
      std::string m_ReadableScheme;

      //! readable scheme
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModelWrapper();

      //! @brief construct from members
      //! @param SP_FUNCTION function to wrap
      //! @param TYPE score type
      //! @param READABLE_SCHEME readable scheme
      ProteinModelWrapper
      (
        const util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, double> > &SP_FUNCTION,
        const Type &TYPE = e_Undefined,
        const std::string &READABLE_SCHEME = "",
        const std::string &SCHEME = ""
      );

      //! @brief Clone function
      //! @return pointer to new ProteinModelWrapper
      ProteinModelWrapper *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get a more readable score scheme
      //! @return a more readable score scheme
      const std::string &GetReadableScheme() const
      {
        return m_ReadableScheme;
      }

      //! @brief get score type
      //! @return score type
      Type GetType() const
      {
        return m_ScoreType;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @return function value of the given argument
      double operator()( const assemble::ProteinModel &ARGUMENT) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write detailed scheme and values to OSTREAM
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::ProteinModel &ARGUMENT,
        std::ostream &OSTREAM
      ) const;

    }; // class ProteinModelWrapper

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_WRAPPER_H_
