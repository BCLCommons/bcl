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
#include "random/bcl_random_histogram_2d_distribution.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace random
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Histogram2DDistribution::s_Instance
    (
      GetObjectInstances().AddInstance( new Histogram2DDistribution())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor which takes a pre-generated histogram2D
    //! @param HISTOGRAM a Histogram2D object which contains the statistics
    Histogram2DDistribution::Histogram2DDistribution( const math::Histogram2D &HISTOGRAM) :
      m_ProbabilityMatrix( HISTOGRAM.GetHistogram())
    {
      m_ProbabilityMatrix.AsVector().SetToSum( 1.0);
    }

    //! @brief virtual copy constructor
    Histogram2DDistribution *Histogram2DDistribution::Clone() const
    {
      return new Histogram2DDistribution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @brief the class name as const &std::string
    const std::string &Histogram2DDistribution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief  determines a biased random case from a probability distribution
    //! @return a randomly drawn 1D position from the probability distribution
    const size_t Histogram2DDistribution::DetermineRandomCase() const
    {
      // initiate a random number from 0 to 1.0
      double temp( GetGlobalRandom().Random< double>( 1.0));

      // initialize the position holder
      size_t pos( 0);

      for
      (
        // initialize the pointers inside the distribution
        const double *ptr( m_ProbabilityMatrix.Begin()), *ptr_end( m_ProbabilityMatrix.End());
        ptr != ptr_end;
        ++ptr, ++pos
      )
      {
        // subtract the random number generated from the value pointed to in m_ProbabilityMatrix
        temp -= *ptr;

        // if the difference is negative or zero, then return the current position in the matrix
        // because the larger the probability the greater the chance of causing this difference to be negative
        if( temp <= 0)
        {
          return pos;
        }
      }
      BCL_Exit( "The probabilities are not set correctly or distribution is empty!!!", -1);
      return pos;
    }

    //! @brief Determines a random 2D object according to the supplied probabilities
    //! @return VectorND< 2, size_t> of the determined random case
    const storage::VectorND< 2, size_t> Histogram2DDistribution::DetermineRandomCase2D() const
    {
      // initialize the 1D case to return a position
      const size_t position( DetermineRandomCase());

      // calculate the 2D coordinates given the 1D position
      const size_t x( position % m_ProbabilityMatrix.GetNumberCols());
      const size_t y( position / m_ProbabilityMatrix.GetNumberRows());

      return storage::VectorND< 2, size_t>( x, y);
    }

    //! @brief Return a random object according to the supplied probabilities
    //! @param X_BOUNDARY the first boundary in the x-direction
    //! @param Y_BOUNDARY the first boundary in the y-direction
    //! @param X_BINSIZE the size of the bin in the x-direction
    //! @param Y_BINSIZE the size of the bin in the y-direction
    //! @return VectorND<2, double> containing the coordinates of the determined random case transformed into the
    //! @return proper units
    const storage::VectorND< 2, double> Histogram2DDistribution::DetermineRandomCase2D
    (
      const double &X_BOUNDARY,
      const double &Y_BOUNDARY,
      const double &X_BINSIZE,
      const double &Y_BINSIZE
    ) const
    {
      // invoke the DetermineRandomCase2D to be used for the conversion
      const storage::VectorND< 2, size_t> position_2d( DetermineRandomCase2D());

      // given the 2D row/column position, calculate the "numerical" position
      const double x( X_BOUNDARY + ( position_2d.First() + 0.5) * X_BINSIZE);
      const double y( Y_BOUNDARY + ( position_2d.Second() + 0.5) * Y_BINSIZE);

      return storage::VectorND< 2, double>( x, y);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read the ProbabilityMatrix from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Histogram2DDistribution::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ProbabilityMatrix, ISTREAM);

      //end
      return ISTREAM;
      }

    //! @brief write ProbabilityMatrix to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &Histogram2DDistribution::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ProbabilityMatrix, OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  } // namespace random
} // namespace bcl
