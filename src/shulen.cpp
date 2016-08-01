/***** shulen.c ***************************************************
 * Description: Compute the null distribution of shustring lengths
 *              in a query/sbjct comparison given the G/C content
 *              of query and sbjct.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:03:03 2004.
 * License: GNU General Public
 *****************************************************************/
#include <cmath>
#include <cfloat>
#include "util.h"

long long bctable[] = {
    6LL, 10LL, 15LL, 20LL, 21LL, 35LL, 28LL, 56LL, 70LL, 36LL, 84LL, 126LL,
    45LL, 120LL, 210LL, 252LL, 55LL, 165LL, 330LL, 462LL, 66LL, 220LL,
    495LL, 792LL, 924LL, 78LL, 286LL, 715LL, 1287LL, 1716LL, 91LL, 364LL,
    1001LL, 2002LL, 3003LL, 3432LL, 105LL, 455LL, 1365LL, 3003LL, 5005LL,
    6435LL, 120LL, 560LL, 1820LL, 4368LL, 8008LL, 11440LL, 12870LL, 136LL,
    680LL, 2380LL, 6188LL, 12376LL, 19448LL, 24310LL, 153LL, 816LL, 3060LL,
    8568LL, 18564LL, 31824LL, 43758LL, 48620LL, 171LL, 969LL, 3876LL,
    11628LL, 27132LL, 50388LL, 75582LL, 92378LL, 190LL, 1140LL, 4845LL,
    15504LL, 38760LL, 77520LL, 125970LL, 167960LL, 184756LL, 210LL, 1330LL,
    5985LL, 20349LL, 54264LL, 116280LL, 203490LL, 293930LL, 352716LL, 231LL,
    1540LL, 7315LL, 26334LL, 74613LL, 170544LL, 319770LL, 497420LL,
    646646LL, 705432LL, 253LL, 1771LL, 8855LL, 33649LL, 100947LL, 245157LL,
    490314LL, 817190LL, 1144066LL, 1352078LL, 276LL, 2024LL, 10626LL,
    42504LL, 134596LL, 346104LL, 735471LL, 1307504LL, 1961256LL, 2496144LL,
    2704156LL, 300LL, 2300LL, 12650LL, 53130LL, 177100LL, 480700LL,
    1081575LL, 2042975LL, 3268760LL, 4457400LL, 5200300LL, 325LL, 2600LL,
    14950LL, 65780LL, 230230LL, 657800LL, 1562275LL, 3124550LL, 5311735LL,
    7726160LL, 9657700LL, 10400600LL, 351LL, 2925LL, 17550LL, 80730LL,
    296010LL, 888030LL, 2220075LL, 4686825LL, 8436285LL, 13037895LL,
    17383860LL, 20058300LL, 378LL, 3276LL, 20475LL, 98280LL, 376740LL,
    1184040LL, 3108105LL, 6906900LL, 13123110LL, 21474180LL, 30421755LL,
    37442160LL, 40116600LL, 406LL, 3654LL, 23751LL, 118755LL, 475020LL,
    1560780LL, 4292145LL, 10015005LL, 20030010LL, 34597290LL, 51895935LL,
    67863915LL, 77558760LL, 435LL, 4060LL, 27405LL, 142506LL, 593775LL,
    2035800LL, 5852925LL, 14307150LL, 30045015LL, 54627300LL, 86493225LL,
    119759850LL, 145422675LL, 155117520LL, 465LL, 4495LL, 31465LL, 169911LL,
    736281LL, 2629575LL, 7888725LL, 20160075LL, 44352165LL, 84672315LL,
    141120525LL, 206253075LL, 265182525LL, 300540195LL, 496LL, 4960LL,
    35960LL, 201376LL, 906192LL, 3365856LL, 10518300LL, 28048800LL,
    64512240LL, 129024480LL, 225792840LL, 347373600LL, 471435600LL,
    565722720LL, 601080390LL, 528LL, 5456LL, 40920LL, 237336LL, 1107568LL,
    4272048LL, 13884156LL, 38567100LL, 92561040LL, 193536720LL, 354817320LL,
    573166440LL, 818809200LL, 1037158320LL, 1166803110LL, 561LL, 5984LL,
    46376LL, 278256LL, 1344904LL, 5379616LL, 18156204LL, 52451256LL,
    131128140LL, 286097760LL, 548354040LL, 927983760LL, 1391975640LL,
    1855967520LL, 2203961430LL, 2333606220LL, 595LL, 6545LL, 52360LL,
    324632LL, 1623160LL, 6724520LL, 23535820LL, 70607460LL, 183579396LL,
    417225900LL, 834451800LL, 1476337800LL, 2319959400LL, 3247943160LL,
    4059928950LL, 4537567650LL, 630LL, 7140LL, 58905LL, 376992LL, 1947792LL,
    8347680LL, 30260340LL, 94143280LL, 254186856LL, 600805296LL,
    1251677700LL, 2310789600LL, 3796297200LL, 5567902560LL, 7307872110LL,
    8597496600LL, 9075135300LL, 666LL, 7770LL, 66045LL, 435897LL, 2324784LL,
    10295472LL, 38608020LL, 124403620LL, 348330136LL, 854992152LL,
    1852482996LL, 3562467300LL, 6107086800LL, 9364199760LL, 12875774670LL,
    15905368710LL, 17672631900LL, 703LL, 8436LL, 73815LL, 501942LL,
    2760681LL, 12620256LL, 48903492LL, 163011640LL, 472733756LL,
    1203322288LL, 2707475148LL, 5414950296LL, 9669554100LL, 15471286560LL,
    22239974430LL, 28781143380LL, 33578000610LL, 35345263800LL, 741LL,
    9139LL, 82251LL, 575757LL, 3262623LL, 15380937LL, 61523748LL,
    211915132LL, 635745396LL, 1676056044LL, 3910797436LL, 8122425444LL,
    15084504396LL, 25140840660LL, 37711260990LL, 51021117810LL,
    62359143990LL, 68923264410LL, 780LL, 9880LL, 91390LL, 658008LL,
    3838380LL, 18643560LL, 76904685LL, 273438880LL, 847660528LL,
    2311801440LL, 5586853480LL, 12033222880LL, 23206929840LL, 40225345056LL,
    62852101650LL, 88732378800LL, 113380261800LL, 131282408400LL,
    137846528820LL, 820LL, 10660LL, 101270LL, 749398LL, 4496388LL,
    22481940LL, 95548245LL, 350343565LL, 1121099408LL, 3159461968LL,
    7898654920LL, 17620076360LL, 35240152720LL, 63432274896LL,
    103077446706LL, 151584480450LL, 202112640600LL, 244662670200LL,
    269128937220LL, 861LL, 11480LL, 111930LL, 850668LL, 5245786LL,
    26978328LL, 118030185LL, 445891810LL, 1471442973LL, 4280561376LL,
    11058116888LL, 25518731280LL, 52860229080LL, 98672427616LL,
    166509721602LL, 254661927156LL, 353697121050LL, 446775310800LL,
    513791607420LL, 538257874440LL, 903LL, 12341LL, 123410LL, 962598LL,
    6096454LL, 32224114LL, 145008513LL, 563921995LL, 1917334783LL,
    5752004349LL, 15338678264LL, 36576848168LL, 78378960360LL,
    151532656696LL, 265182149218LL, 421171648758LL, 608359048206LL,
    800472431850LL, 960566918220LL, 1052049481860LL, 946LL, 13244LL,
    135751LL, 1086008LL, 7059052LL, 38320568LL, 177232627LL, 708930508LL,
    2481256778LL, 7669339132LL, 21090682613LL, 51915526432LL,
    114955808528LL, 229911617056LL, 416714805914LL, 686353797976LL,
    1029530696964LL, 1408831480056LL, 1761039350070LL, 2012616400080LL,
    2104098963720LL, 990LL, 14190LL, 148995LL, 1221759LL, 8145060LL,
    45379620LL, 215553195LL, 886163135LL, 3190187286LL, 10150595910LL,
    28760021745LL, 73006209045LL, 166871334960LL, 344867425584LL,
    646626422970LL, 1103068603890LL, 1715884494940LL, 2438362177020LL,
    3169870830126LL, 3773655750150LL, 4116715363800LL, 1035LL, 15180LL,
    163185LL, 1370754LL, 9366819LL, 53524680LL, 260932815LL, 1101716330LL,
    4076350421LL, 13340783196LL, 38910617655LL, 101766230790LL,
    239877544005LL, 511738760544LL, 991493848554LL, 1749695026860LL,
    2818953098830LL, 4154246671960LL, 5608233007146LL, 6943526580276LL,
    7890371113950LL, 8233430727600LL, 1081LL, 16215LL, 178365LL, 1533939LL,
    10737573LL, 62891499LL, 314457495LL, 1362649145LL, 5178066751LL,
    17417133617LL, 52251400851LL, 140676848445LL, 341643774795LL,
    751616304549LL, 1503232609098LL, 2741188875414LL, 4568648125690LL,
    6973199770790LL, 9762479679106LL, 12551759587422LL, 14833897694226LL,
    16123801841550LL, 1128LL, 17296LL, 194580LL, 1712304LL, 12271512LL,
    73629072LL, 377348994LL, 1677106640LL, 6540715896LL, 22595200368LL,
    69668534468LL, 192928249296LL, 482320623240LL, 1093260079344LL,
    2254848913647LL, 4244421484512LL, 7309837001104LL, 11541847896480LL,
    16735679449896LL, 22314239266528LL, 27385657281648LL, 30957699535776LL,
    32247603683100LL, 1176LL, 18424LL, 211876LL, 1906884LL, 13983816LL,
    85900584LL, 450978066LL, 2054455634LL, 8217822536LL, 29135916264LL,
    92263734836LL, 262596783764LL, 675248872536LL, 1575580702584LL,
    3348108992991LL, 6499270398159LL, 11554258485616LL, 18851684897584LL,
    28277527346376LL, 39049918716424LL, 49699896548176LL, 58343356817424LL,
    63205303218876LL, 1225LL, 19600LL, 230300LL, 2118760LL, 15890700LL,
    99884400LL, 536878650LL, 2505433700LL, 10272278170LL, 37353738800LL,
    121399651100LL, 354860518600LL, 937845656300LL, 2250829575120LL,
    4923689695575LL, 9847379391150LL, 18053528883775LL, 30405943383200LL,
    47129212243960LL, 67327446062800LL, 88749815264600LL, 108043253365600LL,
    121548660036300LL, 126410606437752LL, 1275LL, 20825LL, 249900LL,
    2349060LL, 18009460LL, 115775100LL, 636763050LL, 3042312350LL,
    12777711870LL, 47626016970LL, 158753389900LL, 476260169700LL,
    1292706174900LL, 3188675231420LL, 7174519270695LL, 14771069086725LL,
    27900908274925LL, 48459472266975LL, 77535155627160LL, 114456658306760LL,
    156077261327400LL, 196793068630200LL, 229591913401900LL,
    247959266474052LL, 1326LL, 22100LL, 270725LL, 2598960LL, 20358520LL,
    133784560LL, 752538150LL, 3679075400LL, 15820024220LL, 60403728840LL,
    206379406870LL, 635013559600LL, 1768966344600LL, 4481381406320LL,
    10363194502115LL, 21945588357420LL, 42671977361650LL, 76360380541900LL,
    125994627894135LL, 191991813933920LL, 270533919634160LL,
    352870329957600LL, 426384982032100LL, 477551179875952LL,
    495918532948104LL, 1378LL, 23426LL, 292825LL, 2869685LL, 22957480LL,
    154143080LL, 886322710LL, 4431613550LL, 19499099620LL, 76223753060LL,
    266783135710LL, 841392966470LL, 2403979904200LL, 6250347750920LL,
    14844575908435LL, 32308782859535LL, 64617565719070LL, 119032357903550LL,
    202355008436035LL, 317986441828055LL, 462525733568080LL,
    623404249591760LL, 779255311989700LL, 903936161908052LL,
    973469712824056LL, 1431LL, 24804LL, 316251LL, 3162510LL, 25827165LL,
    177100560LL, 1040465790LL, 5317936260LL, 23930713170LL, 95722852680LL,
    343006888770LL, 1108176102180LL, 3245372870670LL, 8654327655120LL,
    21094923659355LL, 47153358767970LL, 96926348578605LL, 183649923622620LL,
    321387366339585LL, 520341450264090LL, 780512175396135LL,
    1085929983159840LL, 1402659561581460LL, 1683191473897752LL,
    1877405874732108LL, 1946939425648112LL,
}; /* 676 Entries */

//from: http://stackoverflow.com/questions/11032781/fastest-way-to-generate-binomial-coefficients
long long binomial(int n, int k) {
    int i;
    long long b;

    if (0 == k || n == k) return 1LL;
    if (k > n) return 0LL;

    if (k > (n - k)) k = n - k;
    if (1 == k) return (long long)n;

    if (n <= 54 && k <= 54) {
        return bctable[(((n - 3) * (n - 3)) >> 2) + (k - 2)];
    }
    /* Last resort: actually calculate */
    b = 1LL;
    for (i = 1; i <= k; ++i) {
        b *= (n - (k - i));
        if (b < 0) return -1LL; /* Overflow */
        b /= i;
    }
    return b;
}

bool thresholdReached = false;

double sum(double x, double p, double l) {
  double s = 0;
  double k = 0;
  if (!thresholdReached) {
    for (k = 0; k <= x; k++) {
      // double binom = gsl_sf_lnchoose(x, k);
      double binom = log((double)binomial(x, k));
      double pows = pow(2, x) * pow(p, k) * pow(0.5 - p, x - k) *
                    pow(1 - pow(p, k) * pow(0.5 - p, x - k), l);
      s += exp(log(pows) + binom);
      if (s >= 1.0 - DBL_EPSILON) {
        thresholdReached = true;
        s = 1.0;
      }
    }
  } else {
    s = 1.0;
  }
  return s;
}

double expShulen(double gc, double l) {
  double cp;   /* cumulative probability */
  double p;    /* G/C-content of query */
  double d;    /* maximum divergence */
  double prob; /* probability */
  double m;    /* mean shustring length */
  double x;    /* current shustring length */
  double prevP1, curP1;
  thresholdReached = false; // reset flag

  p = gc;
  cp = 0.0;
  x = 0;
  m = 0.0;
  prevP1 = 0.0;
  d = 1.0 - gc * gc;
  while (cp < 1.0 - DBL_EPSILON) {
    x++;
    curP1 = sum(x, p / 2, l); /* exact formula */
    curP1 *= 1.0 - pow(1.0 - d, x);
    prob = curP1 - prevP1; /* exact probability */
    prevP1 = curP1;
    if (prob < 0)
      prob = 0;
    cp += prob;
    m += (double)x * prob;
  }
  return m;
}
