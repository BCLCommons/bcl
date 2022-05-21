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

#ifndef BCL_SSPRED_METHOD_INTERFACE_H_
#define BCL_SSPRED_METHOD_INTERFACE_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sspred_methods.h"
#include "biol/bcl_biol_environment_types.h"
#include "biol/bcl_biol_ss_types.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically
namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MethodInterface
    //! @brief the interface class for all secondary structure methods to be derived from.
    //! @details This class provided the interface for all various secondary structure prediction method classes to be derived
    //! from. It has no data storage, and it provides the interface for functions for getting three state or nine state
    //! predictions and reading predictions for a single amino acid and for a full sequence. These functions should be
    //! overwritten by the derived methods classes since each secondary structure method has their own format which
    //! requires distinct extracting, reading and writing
    //! functions.
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Jun 3, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MethodInterface :
      public util::ObjectInterface
    {
    public:

    //////////
    // data //
    //////////

      //! @brief returns static predictions vector3D set to coil
      //! @return static predictions vector3D set to coil
      static const linal::Vector3D &GetDefaultPredictionVector();

      //! @brief returns static predictions linal::Matrix set to solution coil
      //! @return static predictions linal::Matrix set to solution  oil
      static const linal::Matrix< double> &GetDefaultPredictionMatrix();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new MethodInterface
      virtual MethodInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief get file extension associated with this Method
      //! @return file extension associated with this Method
      virtual const std::string &GetFileExtension() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief get three state, environment independent secondary structure prediction
      //! @return three state, environment independent secondary structure prediction
      virtual linal::Vector3D GetThreeStatePrediction() const = 0;

      //! @brief get three state, environment prediction
      //! @return three state, environment prediction
      virtual linal::Vector3D GetThreeStateTMPrediction() const;

      //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
      //! @return three state secondary structure prediction
      virtual linal::Matrix< double> GetNineStatePrediction() const = 0;

      //! @brief get whether this method determined the secondary structure / membrane environment from the structure
      //! @return true if this method determined the secondary structure / membrane environment from the structure
      virtual bool GetIsDeterminedFromSturcture() const
      {
        return false;
      }

      //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AMINO_ACID amino acid into which sspredictions will be read
      //! @return std::istream which was read from
      virtual std::istream &ReadPredictionsForAA
      (
        std::istream &ISTREAM,
        biol::AABase &AMINO_ACID
      ) const = 0;

      //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which sspredictions will be read
      //! @return std::istream which was read from
      virtual std::istream &ReadPredictionsForAASequence
      (
        std::istream &ISTREAM,
        biol::AASequence &AA_SEQUENCE
      ) const = 0;

      //! @brief find the SSType with highest prediction and returns it
      //! @return SSType with highest prediction
      virtual biol::SSType GetOneStateSSPrediction() const;

      //! @brief find the TMTypes with highest prediction and returns it
      //! @return TMTYpe with highest prediction
      virtual biol::EnvironmentType GetOneStateTMPrediction() const;

      //! @brief find the SSType-TMType pair with highest prediction and returns it
      //! @return SSType-TMType pair with highest prediction
      virtual storage::Pair< biol::SSType, biol::EnvironmentType> GetOneStateSSTMPrediction() const;

      //! @brief converts a three state prediction to a nine state predictions for given TM_TYPE
      //! @param THREE_STATE three state secondary structure prediction for given TM_TYPE
      //! @param TM_TYPE environment type which THREE_STATE predictions are for
      //! @return nine state prediction which has the provided THREE_STATE for given TM_TYPE and 0's for rest
      static linal::Matrix< double> ConvertThreeStateToNineState
      (
        const linal::Vector3D &THREE_STATE,
        const biol::EnvironmentType &TM_TYPE
      );

      //! @brief converts a nine state prediction to a three state predictions
      //! @param NINE_STATE nine state secondary structure prediction
      //! @return three state prediction constructed from a nine state prediction
      static linal::Vector3D ConvertNineStateToThreeState
      (
        const linal::Matrix< double> &NINE_STATE
      );

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write secondary structure predictions in three-state form to given OSTREAM
      //! @param OSTREAM output stream
      //! @return std::ostream which was written to
      std::ostream &WriteThreeStatePredictions( std::ostream &OSTREAM) const;

      //! @brief write secondary structure predictions in nine-state form to given OSTREAM
      //! @param OSTREAM output stream
      //! @return std::ostream which was written to
      std::ostream &WriteNineStatePredictions( std::ostream &OSTREAM) const;

      //! @brief write secondary structure predictions to the provided OSTREAM
      //! this function outputs three state predictions unless overwritten by the method
      //! @param OSTREAM output stream
      //! @return std::ostream which was written to
      virtual std::ostream &WritePredictions( std::ostream &OSTREAM) const;

    protected:

      //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which secondary structure predictions will be read
      //! @param SS_METHOD Method to be read
      //! @return std::istream which was read from
      static std::istream &ReadStandardPredictionsForAASequence
      (
        std::istream &ISTREAM,
        biol::AASequence &AA_SEQUENCE,
        const Method &SS_METHOD
      );

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief creates a default prediction matrix with coil,solution set to 1, everything else to 0
      //! @return a default prediction matrix with coil,solution set to 1, everything else to 0
      static linal::Matrix< double> CreateDefaultPredictionMatrix();

    }; // class MethodInterface

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_METHOD_INTERFACE_H_
