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

#ifndef BCL_OPENCL_DENSITY_FIT_PROTEIN_MINIMIZER_POWELL_H_
#define BCL_OPENCL_DENSITY_FIT_PROTEIN_MINIMIZER_POWELL_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_coordinates_transformer.h"
#include "bcl_opencl_dataset_min_max.h"
#include "bcl_opencl_density_simulate_gaussian_sphere.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_protein_agreement_ccc.h"
#include "bcl_opencl_vector.h"
#include "density/bcl_density_fit_protein_minimizer_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DensityFitProteinMinimzerPowell
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_density_density_fit_protein_minimzer_powell.cpp @endlink
    //! @author woetzen
    //! @date Mar 12, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DensityFitProteinMinimzerPowell :
      public density::FitProteinMinimizerInterface
    {

    private:

    ///////////
    // types //
    ///////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class PositionCorrelation
      //!
      //! @author woetzen
      //! @date Dec 5, 2010
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class PositionCorrelation :
        public math::FunctionInterfaceSerializable< linal::Vector< double>, double>
      {

      private:

      //////////
      // data //
      //////////

        CommandQueue                     m_CommandQueue;
        DensitySimulateGaussianSphere    m_Simulator;
        ProteinAgreementCCC              m_Correlation;
        DataSetMinMax< double>           m_MinMax;
        CoordinateTransformer< double>   m_Transformer;
        Matrix< double>                  m_Atoms;
        size_t                           m_NrAtoms;
        util::SiPtr< const density::Map> m_SpMap;
        Vector< double>                  m_DensityMapBuffer;
        storage::VectorND< 3, size_t>    m_Dimensions;

      public:

        //! single instance of that class
        static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

        //! @brief default constructor
        PositionCorrelation();

        //! @brief Clone function
        //! @return pointer to new PositionCorrelation
        PositionCorrelation *Clone() const;

      /////////////////
      // data access //
      /////////////////

        //! @brief returns class name
        //! @return the class name as const ref std::string
        const std::string &GetClassIdentifier() const;

        //! @brief set the protein model
        //! @param PROTEIN_MODEL
        void SetProtein( const assemble::ProteinModel &PROTEIN_MODEL);

        //! @brief set the map
        //! @param DENSITY_MAP
        void SetDensityMap( const util::SiPtr< const density::Map> &DENSITY_MAP);

        //! @brief set the resolution
        //! @param RESOLUTION
        void SetResolution( const double RESOLUTION);

      ////////////////
      // operations //
      ////////////////

        //! @brief is this class compatible with given command queue
        //! @param COMMAND_QUEUE the command queue this object would operate on
        //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
        bool IsCompatible( const CommandQueue &COMMAND_QUEUE) const;

        //! @brief initialize this class
        //! @brief COMMAND_QUEUE queue to use
        //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
        bool Initialize( const CommandQueue &COMMAND_QUEUE);

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

      }; // class PositionCorrelation

    //////////
    // data //
    //////////

      //! max translation
      double m_MaxTranslation;

      //! max rotation in radians
      double m_MaxRotation;

      //! max number of iterations
      size_t m_MaxIterations;

      //! correlation function used
      mutable util::ShPtr< PositionCorrelation> m_CorrelationFunction;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DensityFitProteinMinimzerPowell();

      //! @brief Clone function
      //! @return pointer to new DensityFitProteinMinimzerPowell
      DensityFitProteinMinimzerPowell *Clone() const;

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
      void SetProteinAgreement( const density::ProteinAgreement &AGREEMENT);

      //! @brief simulator to use
      //! @param DENSITY_SIMULATOR simulator enumerator
      void SetSimulator( const density::Simulator &DENSITY_SIMULATOR);

    ////////////////
    // operations //
    ////////////////

      //! @brief is this class compatible with given command queue
      //! @param COMMAND_QUEUE the command queue this object would operate on
      //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
      bool IsCompatible( const CommandQueue &COMMAND_QUEUE) const;

      //! @brief initialize this class
      //! @brief COMMAND_QUEUE queue to use
      //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
      bool Initialize( const CommandQueue &COMMAND_QUEUE);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator minimizing the position of a protein model within a given density map
      //! @param PROTEIN_MODEL start position of given protein model
      //! @param DENSITY_MAP the density map to fit the PROTEIN_MODEL into
      //! @return the fitted protein model
      assemble::ProteinModel operator()( const assemble::ProteinModel &PROTEIN_MODEL, const density::Map &DENSITY_MAP) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:
  
    }; // class DensityFitProteinMinimzerPowell

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_DENSITY_FIT_PROTEIN_MINIMIZER_POWELL_H_ 
