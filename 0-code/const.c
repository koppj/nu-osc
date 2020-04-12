/***************************************************************************
 * Definition of certain constants                                         *
 ***************************************************************************
 * Author: Joachim Kopp                                                    *
 ***************************************************************************/
#include "nu.h"

const env_param actions[] = 
{
  { "SPECTRUM",               NU_ACTION_SPECTRUM               },
  { "PARAM_SCAN",             NU_ACTION_PARAM_SCAN             },
  { "MCMC",                   NU_ACTION_MCMC                   },
  { "EXPOSURE_SCAN",          NU_ACTION_EXPOSURE_SCAN          },
  { "PROB_TABLE",             NU_ACTION_PROB_TABLE             },
//  { "CHECK_BF",               NU_ACTION_CHECK_BF               },
  { NULL, -1  }
};


const env_param ext_scenarios[] =
{
  { "MB-Thomas",        EXT_MB               },
  { "MB300-Thomas",     EXT_MB_300           },
  { "MBANTI-Thomas",    EXT_MBANTI           },
  { "MBANTI200-Thomas", EXT_MBANTI_200       },
  { "KARMEN",           EXT_KARMEN           },
  { "LSND",             EXT_LSND             },
  { "SBL",              EXT_REACTORS         },
  { "REACTORS",         EXT_REACTORS         },
  { "NOMAD",            EXT_NOMAD            },
  { "CDHS",             EXT_CDHS             },
  { "ATM_TABLE",        EXT_ATM_TABLE        },
  { "ATM_COMP",         EXT_ATM_COMP         },
  { "ATM",              EXT_ATM_COMP         },
  { "DEEPCORE",         EXT_DEEPCORE         },
  { "SOLAR",            EXT_SOLAR            },
  { "MINOS_2016",       EXT_MINOS2016        },
  { "MINOS_2017",       EXT_MINOS2017        },
#ifdef NU_USE_NUSQUIDS
  { "LSND_IVAN",        EXT_LSND_IVAN        },
  { "KARMEN_IVAN",      EXT_KARMEN_IVAN      },
  { "DECAY_KINEMATICS", EXT_DECAY_KINEMATICS },
  { "FREE_STREAMING",   EXT_FREE_STREAMING   },
#endif
  { "MB_JK",            EXT_MB_JK            },
  { NULL, -1 }
};
