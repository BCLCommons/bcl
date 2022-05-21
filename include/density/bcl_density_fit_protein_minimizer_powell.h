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

#ifndef BCL_DENSITY_FIT_PROTEIN_MINIMIZER_POWELL_H_
#define BCL_DENSITY_FIT_PROTEIN_MINIMIZER_POWELL_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_density_fit_protein_minimizer_interface.h"
#include "bcl_density_protein_agreements.h"
#include "bcl_density_simulators.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FitProteinMinimizerPowell
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_namespace_fit_protein_minimizer_powell.cpp @endlink
    //! @author woetzen
    //! @date Mar 12, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FitProteinMinimizerPowell :
      public FitProteinMinimizerInterface
    {

    public:

    ///////////
    // types //
    ///////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class PositionCorrelation
      //! @brief function object that is constructed from a set of atoms and a density map, and returns a correlation for
      //!        a vector containing the relative translation and rotation of the atoms
      //! @details TODO: add an detailed description to this class
      //!
      //! @author woetzen
      //! @date Dec 6, 2010
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class PositionCorrelation :
        public math::FunctionInterfaceSerializable< linal::Vector< double>, double>
      {

      private:

      //////////
      // data //
      //////////

        util::ShPtr< ProteinAgreementInterface>    m_Agreement;
        util::SiPtr< const assemble::ProteinModel> m_Model;

      public:

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief default constructor
        PositionCorrelation
        (
          const util::ShPtr< ProteinAgreementInterface> &AGREEMENT,
          const assemble::ProteinModel &PROTEIN_MODEL
        );

        //! @brief Clone function
        //! @return pointer to new PositionCorrelation
        PositionCorrelation *Clone() const;

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name
        //! @return the class name as const ref std::string
        const std::string &GetClassIdentifier() const;

      ///////////////
      // operators //
      ///////////////

        //! @brief calculate the correlation
        //! @param VECTOR containing 6 elements - three rotations, three translations
        double operator()( const linal::Vector< double> &VECTOR) const;

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

      //////////////////////
      // helper functions //
      //////////////////////

        static util::ShPtr< assemble::ProteinModel> TransformedHardCopy
        (
          const assemble::ProteinModel &PROTEIN_MODEL,
          const linal::Vector< double> &VECTOR
        );

      }; // class PositionCorrelation

    //////////
    // data //
    //////////

      //! resolution
      double m_Resolution;

      //! protein agreement measure to be used
      ProteinAgreement m_ProteinAgreement;

      //! simulator to use
      Simulator m_Simulator;

      //! max translation
      double m_MaxTranslation;

      //! max rotation in radians
      double m_MaxRotation;

      //! max number of iterations
      size_t m_MaxIterations;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FitProteinMinimizerPowell();

      //! @brief Clone function
      //! @return pointer to new FitProteinMinimizerPowell
      FitProteinMinimizerPowell *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief set the resolution of the density map
      //! @param RESOLUTION density map and simulation resolution
      void SetResolution( const double RESOLUTION);

      //! @brief set max translation and rotation
      //! @param MAX_TRANSLATION max translation in any direction for a single iteration
      //! @param MAX_ROTATION max rotation in radians in any direction for a single iteration
      void SetMaxTranslationAndRotation( const double MAX_TRANSLATION, const double MAX_ROTATION);

      //! @brief set the max number of iterations for minimization
      //! @param MAX_NUMBER_ITERATIONS maximum number of iterations for minimization
      void SetMaxIterations( const size_t MAX_NUMBER_ITERATIONS);

      //! @brief set protein agreement measure to be used
      //! @param AGREEMENT protein agreement enumerator
      void SetProteinAgreement( const ProteinAgreement &AGREEMENT);

      //! @brief simulator to use
      //! @param DENSITY_SIMULATOR simulator enumerator
      void SetSimulator( const Simulator &DENSITY_SIMULATOR);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator minimizing the position of a protein model within a given density map
      //! @param PROTEIN_MODEL start position of given protein model
      //! @param DENSITY_MAP the density map to fit the PROTEIN_MODEL into
      //! @return the fitted protein model
      assemble::ProteinModel operator()( const assemble::ProteinModel &PROTEIN_MODEL, const Map &DENSITY_MAP) const;

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

    }; // class FitProteinMinimizerPowell

  } // namespace density
} // namespace bcl

#endif // BCL_DENSITY_FIT_PROTEIN_MINIMIZER_POWELL_H_ 
