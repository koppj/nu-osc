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
  { "EXPOSURE_SCAN",          NU_ACTION_EXPOSURE_SCAN          },
  { "CHECK_BF",               NU_ACTION_CHECK_BF               },
  { NULL, -1  }
};


const env_param ext_scenarios[] =
{
  { "MB",              EXT_MB               },
  { "MBANTI",          EXT_MBANTI           },
  { "KARMEN",          EXT_KARMEN           },
  { "LSND",            EXT_LSND             },
  { "SBL",             EXT_SBL              },
  { "NOMAD",           EXT_NOMAD            },
  { "CDHS",            EXT_CDHS             },
  { "ATM_TABLE",       EXT_ATM_TABLE        },
  { "ATM_COMP",        EXT_ATM_COMP         },
  { NULL, -1 }
};
