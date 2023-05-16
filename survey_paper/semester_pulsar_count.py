blue = ["J2241-5236", "J2145-0750", ',J2222-0137', ',J2248-0101',
    'J2155-3118', 'J2048-1616', 'J2108-3429', 'J2330-2005', 'J2325-0530', 'J0034-0721', 'J2145-0750', 'J2336-01',
    'J2234+2114', 'J2317+2149',
    'J2241-5236',
    'J2330-2005', 'J0152-1637', 'J0206-4028', 'J2354-22', 'J0134-2937', 'J0038-2501',
    'J0030+0451', 'J0034-0721', 'J0034-0534',
    'J2241-5236',
    'J0133-6957', 'J2324-6054',
    'J0206-4028', 'J0255-5304', 'J0051+0423', 'J0152-1637',
]

red = [
    'J0151-0635',
    'J0152-1637',
    'J0255-5304',
    'J0206-4028', 'J0418-4154',
    'J0304+1932',
    'J0450-1248', 'J0452-1759', 'J0459-0210',
    'J0401-7608',
    'J0450-1248', 'J0459-0210',
    'J0452-1759', 'J0520-2553', 'J0450-1248',
    'J0514-4408', 'J0636-4549', 'J0702-4956', 'J0600-5756', 'J0437-4715',
    'J0630-2834', 'J0737-3039A', 'J0636-4549',
    'J0525+1115', 'J0528+2200', 'J0534+2200', 'J0614+2229',
    'J0450-1248', 'J0452-1759', 'J0459-0210', 'J0601-0527',
]

green = [
    'J0624-0424',
    'J0630-2834', 'J0737-3039A', 'J0742-2822', 'J0729-1836',
    'J0742-2822', 'J0749-4247', 'J0820-4114', 'J0820-3921', 'J0835-4510', 'J0837-4135', 'J0838-3947', 'J0924-5302', 'J0942-5552', 'J0955-5304', 'J0959-4809',
    'J0826+2637', 'J0837+0610', 'J0922+0638',
    'J0729-1448', 'J0729-1836', 'J0742-2822', 'J0758-1528', 'J0820-1350', 'J0823+0159',
    'J0837-4135', 'J0856-6137', 'J0905-6019', 'J0924-5302', 'J0924-5814', 'J0942-5657', 'J0842-4851', 'J0907-5157', 'J0902-6325',
    'J0856-6137', 'J0902-6325', 'J0904-7459', 'J0905-6019', 'J0924-5302', 'J0924-5814', 'J0942-5552', 'J0942-5657',
    'J0922+0638', 'J1022+1001', 'J0953+0755',
    'J0908-1739', 'J0820-4114', 'J0742-2822', 'J0820-1350', 'J0835-4510', 'J0837-4135', 'J0855-3331', 'J1012-2337',
    'J0835-4510', 'J0837-4135', 'J0924-5302', 'J0942-5657', 'J0955-5304', 'J0959-4809', 'J1003-4747', 'J1057-5226', 'J1116-4122',
    'J0922+0638', 'J0943+1631', 'J0946+0951', 'J0953+0755',
    'J0908-1739', 'J0944-1354', 'J1018-1642',
    'J1057-5226', 'J1059-5742', 'J1116-4122', 'J1121-5444', 'J1123-4844', 'J1123-6651', 'J1136-5525', 'J1141-6545', 'J1146-6030', 'J1202-5820', 'J1224-6407', 'J1225-5556', 'J1240-4124', 'J1312-5402', 'J1320-5359',
    'J0953+0755', 'J1136+1551',
    'J1012-2337', 'J1018-1642', 'J1034-3224', 'J1041-1942',
]
purple = [
    'J1311-1228', 'J1257-1027',
    'J1313+0931', 'J1300+1240',
    'J1059-5742', 'J1202-5820', 'J1224-6407', 'J1430-6623', 'J1141-6545', 'J1112-6926', 'J1239-6832', 'J1456-6843', 'J1340-6456', 'J1123-6651',
    'J1328-4357', 'J1240-4124', 'J1355-5153', 'J1418-3921', 'J1335-3642', 'J1320-5359',
    'J1311-1228', 'J1332-3032', 'J1418-3921',
    'J1313+0931', 'J1300+1240', 'J1311-1228',
    'J1453-6413', 'J1456-6843', 'J1430-6623', 'J1534-5334', 'J1440-6344', 'J1355-5153', 'J1418-3921',
    'J1543-0620',
    'J1543+0929',
    'J1510-4422', 'J1527-3931', 'J1455-3330', 'J1507-4352', 'J1418-3921', 'J1534-5334', 'J1536-4948',
]

blue = list(set(blue))
red = list(set(red))
green = list(set(green))
purple = list(set(purple))

print(f"Blue:   {len(blue)}")
print(f"Red:    {len(red)} ({len(list(set(blue) & set(red)))})")
print(f"Green:  {len(green)} ({len(list(set(red) & set(green)))})")
print(f"Purple: {len(purple)} ({len(list(set(green) & set(purple)))})")