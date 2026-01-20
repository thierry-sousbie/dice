#ifndef __SOLVER_INTERFACE_PARAMS_HXX__
#define __SOLVER_INTERFACE_PARAMS_HXX__

#include <string>

#include "../internal/namespace.header"

namespace slv
{

	template <class RefinementStatus, class CoarseningStatus>
	class SolverInterfaceParamsT
	{
	public:
		static std::string parserCategory() { return "solver"; }
		static std::string classHeader() { return "solver_interface_params"; }
		static float classVersion() { return 0.15; }
		static float compatibleSinceClassVersion() { return 0.15; }

		double t0;
		double dt;
		double tEnd;
		double timeLimit;
		double timeLimitSafety;
		int noRestart;
		int noRefine;
		int noCoarsen;
		int restartEvery;
		int coarsenEvery;
		int resimulate;
		int resimulateEvery;
		std::string outputDir;
		std::string timingsFileName;
		std::string statisticsFileName;
		std::string stopSignalFileName;
		std::string dumpRestartSignalFileName;

		SolverInterfaceParamsT()
		{
			setDefault();
			serializedVersion = classVersion();
		}

		void setSerializedVersion(float ver)
		{
			serializedVersion = ver;
		}

		void setDefault()
		{
			t0 = 0.0;
			dt = 0.01;
			tEnd = 1.0;
			timeLimit = -1;
			timeLimitSafety = 0.95;
			noRestart = false;
			noRefine = false;
			noCoarsen = true;
			coarsenEvery = 1;
			restartEvery = 50;
			resimulate = 0;
			resimulateEvery = 0;
			outputDir = "";
			timingsFileName = "timings.txt";
			statisticsFileName = "statistics.txt";
			stopSignalFileName = "STOP";
			dumpRestartSignalFileName = "DUMP_RESTART";
		}

		template <class R, class PM>
		void parse(R *reader, PM &paramsManager, bool doNotHandleTime = false)
		{
			if (!doNotHandleTime)
			{
				t0 = paramsManager.get("t0", parserCategory(), t0,
									   reader, PM::FILE_FIRST,
									   "Initial value of time");

				tEnd = paramsManager.get("tEnd", parserCategory(), tEnd,
										 reader, PM::PARSER_FIRST,
										 "Final value of time");

				dt = paramsManager.get("dt", parserCategory(), dt,
									   reader, PM::PARSER_FIRST,
									   "Maximum duration of a time step (can be smaller if adaptive ...)");
			}

			if (RefinementStatus::isEnabled())
			{
				noRefine = paramsManager.get("noRefine", parserCategory(), noRefine,
											 reader, PM::PARSER_FIRST,
											 "Set to prevent mesh refinement");
			}

			if (CoarseningStatus::isEnabled())
			{
				noCoarsen = paramsManager.get("noCoarsen", parserCategory(), noCoarsen,
											  reader, PM::PARSER_FIRST,
											  "Set to prevent mesh coarsening");

				coarsenEvery = paramsManager.get("coarsenEvery", parserCategory(), coarsenEvery,
												 reader, PM::PARSER_FIRST,
												 "Force coarsening every N timesteps");
			}

			restartEvery = paramsManager.get("restartEvery", parserCategory(), restartEvery,
											 reader, PM::PARSER_FIRST,
											 "Write a restart file every N timestep (no restart is written if N<=0)");

			resimulate = paramsManager.get("resimulate", parserCategory(), resimulate,
										   reader, PM::PARSER_FIRST,
										   "Set to start a resimulation after the simulation completes");

			resimulateEvery = paramsManager.get("resimulateEvery", parserCategory(), resimulateEvery,
												reader, PM::PARSER_FIRST,
												"Not implemented");

			outputDir = paramsManager.get("outputDir", parserCategory(), outputDir,
										  reader, PM::PARSER_FIRST,
										  "The directory where all files will be put (default is './', directory must exist)");

			if (outputDir.length() > 0)
				if (outputDir[outputDir.length() - 1] != '/')
					outputDir += "/";

			timingsFileName = paramsManager.get("timingsFileName", parserCategory(), timingsFileName,
												reader, PM::PARSER_FIRST,
												"The name of the file where timings will be stored");

			statisticsFileName = paramsManager.get("statisticsFileName", parserCategory(), statisticsFileName,
												   reader, PM::PARSER_FIRST,
												   "The name of the file where statistics will be stored",
												   serializedVersion > 0.135);

			stopSignalFileName = paramsManager.get("stopSignalFileName", parserCategory(), stopSignalFileName,
												   reader, PM::PARSER_FIRST,
												   "The name of the file to create in the output directory in order to cleanly stop the run.",
												   serializedVersion > 0.115);

			dumpRestartSignalFileName = paramsManager.get("dumpRestartSignalFileName", parserCategory(), dumpRestartSignalFileName,
														  reader, PM::PARSER_FIRST,
														  "The name of the file to create in the output directory in order to force the creatioàn of a restart file on next time step.",
														  serializedVersion > 0.125);

			// Non managed parameters (not saved in restart files)
			typename PM::Parser *parser = paramsManager.getParser();

			timeLimit = parser->get("timeLimit", parserCategory(), timeLimit,
									"CPU time limit (in hours, will stop before spending it). This is an exact limit (<0 for no limit). See also timeLimitSafety");

			timeLimitSafety = parser->get("timeLimitSafety", parserCategory(), timeLimitSafety,
										  "Fraction of the timeLimit while it is safe to keep computing. Leaves (1-timeLimitSafety)*timeLimit for writing the restart file.");

			noRestart = parser->get("noRestart", parserCategory(), noRestart,
									"Set to prevent from dumping any restart file.");
		}

		template <class PP>
		void parse(PP &paramsParser, bool doNotHandleTime = false)
		{
			if (!doNotHandleTime)
			{
				t0 = paramsParser.get("t0", parserCategory(), t0,
									  "Initial value of time");

				tEnd = paramsParser.get("tEnd", parserCategory(), tEnd,
										"Final value of time");

				dt = paramsParser.get("dt", parserCategory(), dt,
									  "Maximum duration of a time step (can be smaller if adaptive ...)");
			}

			if (RefinementStatus::isEnabled())
			{
				noRefine = paramsParser.get("noRefine", parserCategory(), noRefine,
											"Set to prevent mesh refinement");
			}

			if (CoarseningStatus::isEnabled())
			{
				noCoarsen = paramsParser.get("noCoarsen", parserCategory(), noCoarsen,
											 "Set to prevent mesh coarsening");

				coarsenEvery = paramsParser.get("coarsenEvery", parserCategory(), coarsenEvery,
												"Force coarsening every N timesteps");
			}

			restartEvery = paramsParser.get("restartEvery", parserCategory(), restartEvery,
											"Write a restart file every N timestep (no restart is written if N<=0)");

			resimulate = paramsParser.get("resimulate", parserCategory(), resimulate,
										  "Set to start a resimulation after the simulation completes");

			resimulateEvery = paramsParser.get("resimulateEvery", parserCategory(), resimulateEvery,
											   "Not implemented");

			outputDir = paramsParser.get("outputDir", parserCategory(), outputDir,
										 "The directory where all files will be put (default is './', directory must exist)");

			if (outputDir.length() > 0)
				if (outputDir[outputDir.length() - 1] != '/')
					outputDir += "/";

			timingsFileName = paramsParser.get("timingsFileName", parserCategory(), timingsFileName,
											   "The name of the file where timings will be stored");

			statisticsFileName = paramsParser.get("statisticsFileName", parserCategory(), statisticsFileName,
												  "The name of the file where statistics will be stored");

			stopSignalFileName = paramsParser.get("stopSignalFileName", parserCategory(), stopSignalFileName,
												  "The name of the file to create in the output directory in order to force the run to stop at the end of the current time step.");

			dumpRestartSignalFileName = paramsParser.get("dumpRestartSignalFileName", parserCategory(), dumpRestartSignalFileName,
														 "The name of the file to create in the output directory in order to force the creatioàn of a restart file on next time step.");

			timeLimit = paramsParser.get("timeLimit", parserCategory(), timeLimit,
										 "CPU time limit (in hours, will stop before spending it)");

			timeLimitSafety = paramsParser.get("timeLimitSafety", parserCategory(), timeLimitSafety,
											   "Fraction of the timeLimit while it is safe to keep computing.");

			noRestart = paramsParser.get("noRestart", parserCategory(), noRestart,
										 "Prevent from dumping any restart file.");
		}

	private:
		// The version of the class from the file we read from
		float serializedVersion;

	}; // class

} // namespace

#include "../internal/namespace.footer"
#endif
