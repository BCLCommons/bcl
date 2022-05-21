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

#ifndef BCL_MODEL_APPROXIMATOR_KAPPA_NEAREST_NEIGHBOR_H_
#define BCL_MODEL_APPROXIMATOR_KAPPA_NEAREST_NEIGHBOR_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorKappaNearestNeighbor
    //! @brief trains a kappa nearest neighbor (optimizing kappa)
    //!
    //! @see @link example_model_approximator_kappa_nearest_neighbor.cpp @endlink
    //! @author mendenjl, fischea
    //! @date Aug 05, 2011
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorKappaNearestNeighbor :
      public ApproximatorBase
    {

    //////////
    // data //
    //////////

    private:

      //! smallest kappa to test
      size_t m_MinKappa;

      //! kappa currently being tested
      size_t m_Kappa;

      //! largest kappa to test
      size_t m_MaxKappa;

      //! last objective function result
      float m_LastObjectiveFunctionResult;

    public:

      //! static instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorKappaNearestNeighbor();

      //! @brief clone function
      //! @return pointer to new ApproximatorKappaNearestNeighbor
      ApproximatorKappaNearestNeighbor *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief sets the training data set
      //! @param DATA training data set to be set
      void SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA);

      //! @brief returns the current model
      //! @return current model
      util::ShPtr< Interface> GetCurrentModel() const;

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > GetCurrentApproximation() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ApproximatorKappaNearestNeighbor

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_KAPPA_NEAREST_NEIGHBOR_H_
