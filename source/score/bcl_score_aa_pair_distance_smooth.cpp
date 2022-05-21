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
#include "score/bcl_score_aa_pair_distance_smooth.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &AAPairDistanceSmooth::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "aa_distances.histograms");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &AAPairDistanceSmooth::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "aasmooth");

      // end
      return s_default_scheme;

    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param SCHEME scheme to be used
    AAPairDistanceSmooth::AAPairDistanceSmooth
    (
      const std::string &HISTOGRAM_FILENAME, // = GetDefaultHistogramFilename(),
      const std::string &SCHEME// = GetDefaultScheme()
    ) :
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_Scheme( SCHEME),
      m_EnergyFunctionMap(),
      m_DistanceCutoff( 0.0)
    {
      // read the histogram file and store the energy functions
      ReadEnergyFunctionMap();
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAPairDistanceSmooth object that is copied from this one
    AAPairDistanceSmooth *AAPairDistanceSmooth::Clone() const
    {
      return new AAPairDistanceSmooth( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAPairDistanceSmooth::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AAPairDistanceSmooth::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &AAPairDistanceSmooth::GetAlias() const
    {
      static const std::string s_name( "AAPairDistanceSmooth");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairDistanceSmooth::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B
    ) const
    {
      // calculate the distance
      const double distance( biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B));

      // return the score
      return operator()( AMINO_ACID_A, AMINO_ACID_B, distance);
    }

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param DISTANCE distance between the amino acid pair
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairDistanceSmooth::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const double DISTANCE
    ) const
    {
      // check that distance is defined
      if( !util::IsDefined( DISTANCE))
      {
        return 0.0;
      }

      // construct biol::AAType pair
      storage::Pair< biol::AAType, biol::AAType> type_pair( AMINO_ACID_A.GetType(), AMINO_ACID_B.GetType());

      // search the map for this aa type pair and store the iterator
      auto itr_find( m_EnergyFunctionMap.Find( type_pair));

      // if the itr is not valid
      if( itr_find == m_EnergyFunctionMap.End())
      {
        // return undefined
        return double( 0);
      }

      // now call the scoring function for the found energy function
      return itr_find->second->operator()( DISTANCE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAPairDistanceSmooth::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_HistogramFileName, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // read the histogram file and store the energy functions
      ReadEnergyFunctionMap();

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AAPairDistanceSmooth::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_HistogramFileName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairDistanceSmooth::WriteDetailedSchemeAndValues
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      std::ostream &OSTREAM
    ) const
    {
      //write Scheme
      OSTREAM << AMINO_ACID_A.GetSeqID() << '\t'
              << AMINO_ACID_A.GetType()->GetThreeLetterCode() << '\t'
              << AMINO_ACID_B.GetSeqID() << '\t'
              << AMINO_ACID_B.GetType()->GetThreeLetterCode() << '\t'
              << biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B) << '\t'
              << operator()( AMINO_ACID_A, AMINO_ACID_B) << '\n';

      //end
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPairDistanceSmooth::GetSerializer() const
    {
      io::Serializer serializer;

      serializer.SetClassDescription( "Scoring AA pair distances with a fitted function.");
      serializer.AddInitializer
      (
        "histogram filename",
        "path to file where the statistics and in consequence the energy potentials are read from",
        io::Serialization::GetAgent( &m_HistogramFileName),
        GetDefaultHistogramFilename()
      );

      return serializer;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read map of amino acid pair energies based on distance from histogram files
    void AAPairDistanceSmooth::ReadEnergyFunctionMap()
    {
      // reset distance cutoff to 0
      m_DistanceCutoff = 0.0;

      // read file with all histograms for each pair of sstypes
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      std::map< std::pair< biol::AAType, biol::AAType>, util::ShPtr< AAPairDistanceFittedFunction> > energymap;

      // read the two first amino acid types
      std::string tmp_a, tmp_b;
      while( read >> tmp_a >> tmp_b && !read.eof())
      {
        // initialize the variable that are going to be read
        math::Histogram current_aa_distance_histogram;
        storage::Pair< biol::AAType, biol::AAType> aa_type_pair
        (
          biol::GetAATypes().AATypeFromOneLetterCode( tmp_a[0]),
          biol::GetAATypes().AATypeFromOneLetterCode( tmp_b[0])
        );

        // abort if any of the aatypes if unknown
        if
        (
          aa_type_pair.First() == biol::GetAATypes().e_Undefined ||
          aa_type_pair.Second() == biol::GetAATypes().e_Undefined
        )
        {
          // alert user and break
          BCL_MessageCrt
          (
            "undefined AAType found in the histogram " + util::Format()( aa_type_pair)
          );
          break;
        }

        // read the histogram
        read >> current_aa_distance_histogram;

        //create spline for the current distribution and store it in map also as swapped pair
        util::ShPtr< AAPairDistanceFittedFunction> current_energy_function
        (
          new AAPairDistanceFittedFunction( current_aa_distance_histogram)
        );

        m_DistanceCutoff = std::max( m_DistanceCutoff, current_energy_function->GetDistanceCutoff());

        // store the spline in the map
        m_EnergyFunctionMap[ aa_type_pair] = current_energy_function;

        // reverse the pair and insert it
        storage::Pair< biol::AAType, biol::AAType> swapped_aa_type_pair( aa_type_pair.Second(), aa_type_pair.First());
        m_EnergyFunctionMap[ swapped_aa_type_pair] = current_energy_function;
      }

      // close the stream
      io::File::CloseClearFStream( read);
    }

  } // namespace score
} // namespace bcl
