# lumi masks applied when running on data
# produce list with ConvertJSON.py in test/tools

#used by crab
#JSONFILE = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"
#JSONFILE = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"
#JSONFILE = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"
JSONFILE = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt"

#used in local productions
LUMIMASK = cms.untracked.VLuminosityBlockRange( *(
    '315257:1-315257:88',
    '315257:91-315257:92',
    '315259:1-315259:172',
    '315264:32-315264:261',
    '315265:4-315265:58',
    '315267:1-315267:244',
    '315270:1-315270:633',
    '315322:23-315322:118',
    '315322:122-315322:1354',
    '315339:37-315339:654',
    '315357:44-315357:732',
    '315357:736-315357:770',
    '315357:780-315357:831',
    '315361:40-315361:619',
    '315363:1-315363:35',
    '315363:37-315363:47',
    '315363:49-315363:67',
    '315363:69-315363:80',
    '315363:82-315363:90',
    '315366:10-315366:61',
    '315366:67-315366:750',
    '315420:28-315420:920',
    '315420:924-315420:942',
    '315420:954-315420:1748',
    '315488:42-315488:843',
    '315489:1-315489:653',
    '315489:672-315489:709',
    '315490:1-315490:24',
    '315506:13-315506:100',
    '315510:1-315510:345',
    '315512:1-315512:1122',
    '315543:55-315543:171',
    '315555:22-315555:97',
    '315556:1-315556:26',
    '315557:1-315557:279',
    '315640:46-315640:87',
    '315641:1-315641:4',
    '315642:1-315642:92',
    '315644:1-315644:184',
    '315645:1-315645:40',
    '315645:47-315645:390',
    '315645:395-315645:565',
    '315645:567-315645:594',
    '315646:1-315646:1033',
    '315647:1-315647:58',
    '315648:1-315648:110',
    '315689:24-315689:1127',
    '315689:1180-315689:1186',
    '315690:10-315690:654',
    '315702:38-315702:113',
    '315703:1-315703:545',
    '315704:1-315704:61',
    '315705:1-315705:700',
    '315713:35-315713:359',
    '315713:374-315713:385',
    '315713:400-315713:1123',
    '315721:33-315721:50',
    '315721:56-315721:626',
    '315741:34-315741:92',
    '315764:37-315764:309',
    '315770:39-315770:332',
    '315784:29-315784:33',
    '315784:40-315784:156',
    '315784:158-315784:161',
    '315785:1-315785:198',
    '315785:201-315785:305',
    '315786:1-315786:72',
    '315790:1-315790:716',
    '315790:718-315790:922',
    '315800:41-315800:621',
    '315801:1-315801:344',
    '315840:33-315840:1154',
    '315973:39-315973:240',
    '315973:262-315973:914',
    '315974:1-315974:71',
    '316058:42-316058:405',
    '316059:1-316059:321',
    '316059:323-316059:567',
    '316060:1-316060:935',
    '316061:1-316061:23',
    '316061:194-316061:206',
    '316062:1-316062:4',
    '316082:37-316082:407',
    '316110:1-316110:210',
    '316111:1-316111:48',
    '316113:1-316113:64',
    '316114:1-316114:777',
    '316114:779-316114:1562',
    '316153:1-316153:770',
    '316186:38-316186:81',
    '316187:1-316187:1091',
    '316187:1093-316187:1100',
    '316187:1207-316187:2077',
    '316199:33-316199:1197',
    '316200:1-316200:10',
    '316201:1-316201:498',
    '316202:1-316202:403',
    '316216:25-316216:466',
    '316217:1-316217:264',
    '316218:1-316218:1008',
    '316219:1-316219:283',
    '316239:38-316239:626',
    '316240:1-316240:1224',
    '316241:1-316241:325',
    '316271:36-316271:121',
    '316361:22-316361:124',
    '316361:126-316361:131',
    '316361:133-316361:135',
    '316361:137-316361:137',
    '316361:139-316361:142',
    '316361:144-316361:145',
    '316361:147-316361:147',
    '316361:149-316361:159',
    '316361:161-316361:174',
    '316361:176-316361:178',
    '316361:180-316361:189',
    '316361:191-316361:197',
    '316361:199-316361:208',
    '316361:210-316361:223',
    '316362:1-316362:208',
    '316362:210-316362:212',
    '316362:214-316362:225',
    '316362:227-316362:242',
    '316362:244-316362:269',
    '316362:271-316362:319',
    '316362:332-316362:392',
    '316362:394-316362:395',
    '316362:397-316362:402',
    '316362:404-316362:404',
    '316362:406-316362:410',
    '316362:412-316362:412',
    '316362:414-316362:418',
    '316362:420-316362:428',
    '316362:430-316362:450',
    '316363:1-316363:39',
    '316363:41-316363:49',
    '316377:19-316377:19',
    '316377:21-316377:40',
    '316378:1-316378:29',
    '316379:1-316379:70',
    '316380:1-316380:708',
    '316380:714-316380:1213',
    '316455:36-316455:71',
    '316457:1-316457:1454',
    '316469:17-316469:444',
    '316470:1-316470:476',
    '316472:1-316472:70',
    '316472:76-316472:333',
    '316505:44-316505:205',
    '316505:207-316505:921',
    '316505:923-316505:1364',
    '316569:20-316569:703',
    '316569:742-316569:1945',
    '316590:17-316590:526',
    '316613:49-316613:241',
    '316615:1-316615:338',
    '316666:1-316666:981',
    '316667:1-316667:197',
    '316700:46-316700:346',
    '316700:388-316700:397',
    '316701:1-316701:479',
    '316702:1-316702:388',
    '316715:33-316715:45',
    '316716:1-316716:181',
    '316717:1-316717:192',
    '316718:1-316718:311',
    '316719:1-316719:91',
    '316719:100-316719:144',
    '316720:1-316720:182',
    '316721:1-316721:15',
    '316722:1-316722:751',
    '316723:1-316723:64',
    '316758:11-316758:1609',
    '316766:51-316766:1920',
    '316766:1922-316766:2199',
    '316876:34-316876:38',
    '316876:40-316876:644',
    '316877:1-316877:164',
    '316877:171-316877:401',
    '316879:1-316879:156',
    '316928:40-316928:188',
    '316985:33-316985:503',
    '316993:44-316993:254',
    '316994:1-316994:14',
    '316995:1-316995:623',
    '317080:41-317080:66',
    '317087:43-317087:177',
    '317087:213-317087:222',
    '317087:257-317087:852',
    '317089:1-317089:1003',
    '317182:47-317182:63',
    '317182:65-317182:1424',
    '317212:36-317212:175',
    '317213:1-317213:375',
    '317279:43-317279:508',
    '317291:34-317291:824',
    '317292:1-317292:330',
    '317297:1-317297:283',
    '317297:347-317297:760',
    '317319:44-317319:182',
    '317320:1-317320:326',
    '317320:333-317320:411',
    '317320:413-317320:1827',
    '317338:66-317338:107',
    '317339:1-317339:163',
    '317340:1-317340:418',
    '317382:58-317382:128',
    '317383:1-317383:58',
    '317391:39-317391:46',
    '317392:1-317392:1116',
    '317392:1119-317392:1900',
    '317435:1-317435:1397',
    '317438:1-317438:68',
    '317438:71-317438:309',
    '317475:33-317475:89',
    '317475:105-317475:115',
    '317478:1-317478:23',
    '317484:1-317484:448',
    '317484:467-317484:514',
    '317484:519-317484:545',
    '317488:1-317488:844',
    '317527:41-317527:1487',
    '317591:43-317591:334',
    '317626:40-317626:2045',
    '317640:29-317640:829',
    '317641:1-317641:1390',
    '317648:45-317648:139',
    '317649:1-317649:621',
    '317650:1-317650:1304',
    '317661:35-317661:1256',
    '317663:1-317663:858',
    '317683:83-317683:402',
    '317696:38-317696:682',
    '318733:1-318733:33',
    '318828:54-318828:123',
    '318872:16-318872:287',
    '318874:1-318874:320',
    '318876:1-318876:161',
    '318877:1-318877:615',
    '319077:52-319077:92',
    '319337:48-319337:2240',
    '319347:40-319347:690',
    '319348:1-319348:37',
    '319349:1-319349:148',
    '319449:35-319449:559',
    '319449:562-319449:734',
    '319450:1-319450:287',
    '319450:290-319450:683',
    '319456:138-319456:346',
    '319459:1-319459:78',
    '319486:38-319486:103',
    '319503:1-319503:317',
    '319524:36-319524:1459',
    '319526:1-319526:282',
    '319528:1-319528:259',
    '319579:41-319579:3168',
    '319625:17-319625:206',
    '319639:31-319639:1509',
    '319656:51-319656:310',
    '319657:1-319657:167',
    '319658:1-319658:225',
    '319659:1-319659:87',
    '319678:36-319678:294',
    '319687:46-319687:90',
    '319697:47-319697:482',
    '319697:490-319697:490',
    '319698:1-319698:312',
    '319756:44-319756:1966',
    '319840:41-319840:388',
    '319841:1-319841:167',
    '319847:49-319847:51',
    '319848:1-319848:53',
    '319849:1-319849:492',
    '319851:1-319851:4',
    '319853:1-319853:40',
    '319853:47-319853:262',
    '319854:1-319854:225',
    '319908:1-319908:40',
    '319908:43-319908:53',
    '319909:1-319909:7',
    '319910:1-319910:983',
    '319912:1-319912:59',
    '319913:1-319913:56',
    '319914:1-319914:32',
    '319915:1-319915:416',
    '319941:43-319941:298',
    '319942:1-319942:50',
    '319950:38-319950:205',
    '319991:46-319991:882',
    '319992:1-319992:264',
    '319993:1-319993:955',
    '320002:52-320002:192',
    '320006:1-320006:34',
    '320006:36-320006:341',
    '320010:1-320010:330',
    '320011:1-320011:302',
    '320012:1-320012:99',
    '320023:17-320023:292',
    '320024:1-320024:410',
    '320025:1-320025:113',
    '320026:1-320026:204',
    '320038:43-320038:663',
    '320039:1-320039:30',
    '320040:1-320040:737',
    '320059:1-320059:105',
    '320060:1-320060:42',
    '320061:1-320061:49',
    '320062:1-320062:21',
    '320063:1-320063:64',
    '320064:1-320064:200',
    '320065:1-320065:920',
    '320673:35-320673:901',
    '320674:1-320674:599',
    '320688:49-320688:531',
    '320712:39-320712:242',
    '320757:51-320757:382',
    '320804:46-320804:1274',
    '320807:1-320807:7',
    '320809:1-320809:716',
    '320821:41-320821:221',
    '320822:1-320822:523',
    '320823:1-320823:360',
    '320824:1-320824:1051',
    '320838:93-320838:357',
    '320840:1-320840:471',
    '320841:1-320841:205',
    '320853:41-320853:369',
    '320854:1-320854:125',
    '320855:1-320855:565',
    '320856:1-320856:159',
    '320857:1-320857:272',
    '320858:1-320858:230',
    '320859:1-320859:40',
    '320887:49-320887:321',
    '320888:1-320888:26',
    '320916:2-320916:25',
    '320917:1-320917:1926',
    '320920:1-320920:178',
    '320933:40-320933:214',
    '320934:1-320934:831',
    '320936:1-320936:407',
    '320941:1-320941:93',
    '320980:44-320980:142',
    '320995:26-320995:214',
    '320996:1-320996:380',
    '321004:39-321004:188',
    '321005:1-321005:61',
    '321006:1-321006:162',
    '321007:1-321007:831',
    '321009:1-321009:85',
    '321010:1-321010:342',
    '321011:1-321011:213',
    '321012:1-321012:35',
    '321012:190-321012:201',
    '321051:58-321051:1179',
    '321055:1-321055:302',
    '321055:304-321055:326',
    '321055:328-321055:340',
    '321055:368-321055:759',
    '321067:39-321067:225',
    '321067:232-321067:639',
    '321068:1-321068:715',
    '321069:1-321069:313',
    '321119:45-321119:214',
    '321121:1-321121:47',
    '321122:1-321122:395',
    '321124:1-321124:819',
    '321126:1-321126:493',
    '321134:33-321134:70',
    '321138:1-321138:741',
    '321140:1-321140:798',
    '321149:35-321149:1424',
    '321149:1426-321149:1476',
    '321149:1478-321149:1553',
    '321149:1558-321149:1576',
    '321149:1578-321149:1588',
    '321149:1591-321149:1743',
    '321165:1-321165:8',
    '321166:1-321166:10',
    '321167:1-321167:141',
    '321167:143-321167:143',
    '321167:145-321167:510',
    '321167:512-321167:552',
    '321167:554-321167:691',
    '321167:693-321167:923',
    '321177:38-321177:74',
    '321177:77-321177:214',
    '321177:216-321177:232',
    '321177:234-321177:247',
    '321177:249-321177:321',
    '321177:323-321177:365',
    '321177:367-321177:455',
    '321178:5-321178:78',
    '321218:49-321218:962',
    '321219:1-321219:934',
    '321221:1-321221:40',
    '321230:41-321230:124',
    '321231:1-321231:59',
    '321232:1-321232:30',
    '321233:1-321233:727',
    '321262:1-321262:4',
    '321283:48-321283:357',
    '321294:1-321294:62',
    '321295:1-321295:307',
    '321295:309-321295:316',
    '321295:318-321295:384',
    '321295:390-321295:394',
    '321295:396-321295:604',
    '321295:606-321295:616',
    '321295:619-321295:646',
    '321295:649-321295:690',
    '321295:693-321295:754',
    '321296:1-321296:24',
    '321296:34-321296:41',
    '321296:44-321296:67',
    '321305:20-321305:2600',
    '321305:2605-321305:2651',
    '321311:1-321311:10',
    '321312:1-321312:768',
    '321313:1-321313:408',
    '321393:1-321393:127',
    '321393:134-321393:148',
    '321396:1-321396:1475',
    '321397:1-321397:365',
    '321414:31-321414:1283',
    '321415:1-321415:804',
    '321431:30-321431:189',
    '321432:1-321432:47',
    '321433:1-321433:125',
    '321434:1-321434:642',
    '321436:1-321436:710',
    '321457:43-321457:451',
    '321457:453-321457:1888',
    '321461:1-321461:149',
    '321475:50-321475:518',
    '321475:526-321475:2084',
    '321710:1-321710:57',
    '321712:1-321712:2',
    '321712:16-321712:54',
    '321712:57-321712:115',
    '321712:117-321712:263',
    '321730:2-321730:257',
    '321730:259-321730:291',
    '321732:1-321732:127',
    '321732:129-321732:181',
    '321732:185-321732:189',
    '321732:192-321732:245',
    '321732:248-321732:252',
    '321732:254-321732:373',
    '321732:375-321732:381',
    '321732:386-321732:386',
    '321732:389-321732:392',
    '321732:395-321732:424',
    '321732:426-321732:432',
    '321732:434-321732:448',
    '321732:450-321732:452',
    '321732:454-321732:459',
    '321732:467-321732:586',
    '321732:589-321732:680',
    '321732:682-321732:686',
    '321732:689-321732:903',
    '321732:905-321732:973',
    '321732:975-321732:1448',
    '321735:1-321735:146',
    '321755:33-321755:361',
    '321755:363-321755:470',
    '321755:472-321755:473',
    '321755:475-321755:487',
    '321755:489-321755:729',
    '321758:1-321758:47',
    '321758:49-321758:75',
    '321758:77-321758:121',
    '321758:128-321758:130',
    '321758:146-321758:148',
    '321758:151-321758:155',
    '321758:161-321758:165',
    '321758:168-321758:189',
    '321760:1-321760:171',
    '321760:175-321760:205',
    '321760:207-321760:238',
    '321760:240-321760:258',
    '321760:260-321760:420',
    '321760:422-321760:520',
    '321760:526-321760:586',
    '321760:588-321760:593',
    '321760:598-321760:602',
    '321760:604-321760:607',
    '321760:613-321760:716',
    '321760:719-321760:721',
    '321760:727-321760:788',
    '321760:794-321760:818',
    '321760:822-321760:824',
    '321760:828-321760:830',
    '321760:834-321760:836',
    '321760:840-321760:841',
    '321760:845-321760:855',
    '321773:11-321773:14',
    '321773:25-321773:35',
    '321773:39-321773:52',
    '321773:54-321773:79',
    '321774:1-321774:12',
    '321774:14-321774:52',
    '321774:54-321774:119',
    '321775:1-321775:12',
    '321775:14-321775:14',
    '321776:1-321776:12',
    '321776:15-321776:19',
    '321776:30-321776:45',
    '321777:1-321777:81',
    '321777:83-321777:169',
    '321777:174-321777:176',
    '321777:192-321777:207',
    '321778:8-321778:150',
    '321780:1-321780:332',
    '321780:336-321780:338',
    '321780:342-321780:346',
    '321780:351-321780:357',
    '321780:359-321780:360',
    '321780:362-321780:371',
    '321780:374-321780:383',
    '321780:392-321780:412',
    '321780:414-321780:420',
    '321780:422-321780:493',
    '321780:496-321780:499',
    '321780:502-321780:503',
    '321780:505-321780:508',
    '321780:517-321780:518',
    '321781:6-321781:37',
    '321781:53-321781:56',
    '321781:58-321781:66',
    '321781:69-321781:69',
    '321781:77-321781:180',
    '321781:186-321781:209',
    '321781:212-321781:265',
    '321781:269-321781:274',
    '321781:276-321781:290',
    '321781:293-321781:312',
    '321781:316-321781:410',
    '321781:412-321781:427',
    '321813:32-321813:352',
    '321815:1-321815:23',
    '321817:1-321817:536',
    '321818:1-321818:690',
    '321820:1-321820:214',
    '321831:25-321831:781',
    '321832:1-321832:389',
    '321832:403-321832:510',
    '321833:1-321833:407',
    '321834:1-321834:333',
    '321879:39-321879:47',
    '321879:50-321879:52',
    '321879:55-321879:68',
    '321879:71-321879:73',
    '321879:77-321879:89',
    '321879:93-321879:95',
    '321879:99-321879:111',
    '321879:114-321879:116',
    '321879:120-321879:132',
    '321879:136-321879:138',
    '321879:141-321879:154',
    '321879:157-321879:159',
    '321879:163-321879:175',
    '321879:178-321879:181',
    '321879:185-321879:197',
    '321879:200-321879:202',
    '321879:207-321879:218',
    '321879:222-321879:356',
    '321880:1-321880:41',
    '321880:44-321880:132',
    '321887:54-321887:948',
    '321908:43-321908:472',
    '321909:1-321909:208',
    '321909:210-321909:1654',
    '321917:4-321917:156',
    '321917:164-321917:808',
    '321919:1-321919:6',
    '321933:43-321933:232',
    '321933:235-321933:326',
    '321960:18-321960:47',
    '321961:1-321961:354',
    '321973:37-321973:746',
    '321973:748-321973:968',
    '321973:972-321973:1253',
    '321975:1-321975:866',
    '321988:45-321988:996',
    '321988:1106-321988:1486',
    '321990:1-321990:471',
    '322013:14-322013:22',
    '322014:1-322014:17',
    '322022:42-322022:185',
    '322022:201-322022:1805',
    '322040:32-322040:70',
    '322057:38-322057:58',
    '322068:51-322068:724',
    '322079:39-322079:200',
    '322079:216-322079:393',
    '322079:409-322079:428',
    '322106:48-322106:871',
    '322113:48-322113:159',
    '322118:1-322118:516',
    '322118:530-322118:874',
    '322179:43-322179:820',
    '322179:823-322179:1783',
    '322201:39-322201:266',
    '322204:1-322204:280',
    '322204:282-322204:301',
    '322204:303-322204:331',
    '322204:337-322204:1143',
    '322222:1-322222:526',
    '322252:42-322252:1586',
    '322317:48-322317:101',
    '322319:1-322319:163',
    '322322:1-322322:170',
    '322322:267-322322:1205',
    '322324:1-322324:416',
    '322332:37-322332:1055',
    '322348:40-322348:1505',
    '322355:36-322355:137',
    '322356:1-322356:779',
    '322381:45-322381:577',
    '322407:46-322407:582',
    '322430:46-322430:501',
    '322431:59-322431:1166',
    '322480:60-322480:408',
    '322492:1-322492:1386',
    '322510:37-322510:45',
    '322599:43-322599:294',
    '322602:1-322602:69',
    '322602:72-322602:72',
    '322603:1-322603:10',
    '322605:1-322605:280',
    '322617:1-322617:601',
    '322625:41-322625:484',
    '322625:492-322625:1167',
    '322633:1-322633:249',
    '323414:1-323414:46',
    '323423:1-323423:136',
    '323470:38-323470:172',
    '323470:176-323470:218',
    '323470:223-323470:266',
    '323471:1-323471:238',
    '323472:1-323472:64',
    '323473:1-323473:227',
    '323474:1-323474:355',
    '323475:1-323475:77',
    '323487:42-323487:177',
    '323487:184-323487:498',
    '323488:1-323488:514',
    '323488:555-323488:734',
    '323488:738-323488:793',
    '323492:1-323492:33',
    '323493:1-323493:144',
    '323495:1-323495:187',
    '323524:25-323524:561',
    '323525:1-323525:91',
    '323525:97-323525:1126',
    '323526:1-323526:248',
    '323526:253-323526:466',
    '323693:38-323693:151',
    '323696:1-323696:257',
    '323702:1-323702:808',
    '323725:18-323725:346',
    '323726:1-323726:60',
    '323727:1-323727:83',
    '323727:88-323727:677',
    '323727:682-323727:813',
    '323727:819-323727:822',
    '323727:826-323727:987',
    '323755:27-323755:815',
    '323755:818-323755:823',
    '323755:826-323755:826',
    '323755:828-323755:830',
    '323755:833-323755:861',
    '323755:864-323755:964',
    '323775:38-323775:81',
    '323775:84-323775:171',
    '323778:1-323778:934',
    '323790:45-323790:948',
    '323794:1-323794:68',
    '323841:46-323841:510',
    '323857:1-323857:357',
    '323940:49-323940:1567',
    '323954:1-323954:77',
    '323976:31-323976:85',
    '323978:1-323978:73',
    '323980:1-323980:202',
    '323983:1-323983:188',
    '323997:1-323997:498',
    '324021:44-324021:819',
    '324022:1-324022:554',
    '324077:54-324077:710',
    '324077:712-324077:753',
    '324201:20-324201:834',
    '324201:837-324201:1385',
    '324202:1-324202:240',
    '324205:1-324205:163',
    '324206:1-324206:149',
    '324207:1-324207:34',
    '324209:1-324209:142',
    '324237:33-324237:236',
    '324245:23-324245:1681',
    '324293:39-324293:1440',
    '324293:1442-324293:2176',
    '324293:2178-324293:2342',
    '324315:1-324315:200',
    '324315:203-324315:204',
    '324318:1-324318:332',
    '324420:1-324420:625',
    '324729:1-324729:193',
    '324747:63-324747:1139',
    '324764:1-324764:150',
    '324765:1-324765:481',
    '324769:1-324769:328',
    '324772:1-324772:165',
    '324785:77-324785:664',
    '324791:1-324791:1217',
    '324835:40-324835:230',
    '324835:302-324835:369',
    '324840:1-324840:96',
    '324841:1-324841:1347',
    '324846:1-324846:151',
    '324846:154-324846:517',
    '324878:62-324878:111',
    '324878:113-324878:175',
    '324878:180-324878:1800',
    '324897:30-324897:170',
    '324970:1-324970:425',
    '324970:428-324970:598',
    '324970:606-324970:632',
    '324970:634-324970:1529',
    '324970:1532-324970:2195',
    '324980:39-324980:917',
    '324980:919-324980:954',
    '324980:956-324980:968',
    '324980:1005-324980:1042',
    '324980:1044-324980:2340',
    '324997:29-324997:150',
    '324998:1-324998:368',
    '324999:1-324999:14',
    '325000:1-325000:371',
    '325001:1-325001:105',
    '325001:108-325001:171',
    '325001:173-325001:595',
    '325022:45-325022:1594',
    '325057:42-325057:383',
    '325097:40-325097:96',
    '325098:1-325098:8',
    '325099:1-325099:394',
    '325100:1-325100:254',
    '325101:1-325101:462',
    '325101:464-325101:485',
    '325110:1-325110:21',
    '325117:1-325117:533',
    '325159:48-325159:266',
    '325168:1-325168:21',
    '325169:1-325169:23',
    '325170:1-325170:692',
    '325170:694-325170:1205',
    '325172:1-325172:267',
    '325172:269-325172:485',
))
