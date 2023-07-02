import re

airsport_areas = {
    'Oslo TMA' : ['STARMOEN', 'HOKKSUND', 'EGGEMOEN', 'EINA', 'SUNNVOLLEN', 'HVITTINGFOSS'],
    'Farris TMA' : ['HVITTINGFOSS', 'GVARV', 'LUNDE', 'DRANGEDAL', 'BOE', 'BØ', 'TOKKE'],
    'WestCoast TMA' : ['KVAM TRANSITAREA', 'GULLFJELLET'],
    'Værnes TMA' : ['NIDAROS', 'GAULDAL', 'MERÅKER'],
#    'Polaris CTA' : ['Bjorli Wave', 'Lesja Wave'],
    }

airspace_frequencies = {
    'Airwork A' : '118.475',
    'Airwork B' : '118.475',
    'Airwork C2' : '118.475',
    'Airwork E' : '118.475',
    'Airwork F' : '118.475',    
    'Alta TMA' : '120.40',
    'Anda TIZ' : '119.80',
    'Andøya TMA' : '118.20',
    'Banak CTR' : '118.90',
    'Bardufoss TMA' : '118.80',
    'Bergen Lufthavn Flesland' : '119.10',
    'Berlevåg TIZ' : '120.10',
    'Bjorli Wave 125.70' : '125.70',
    'Bodø CTR' : '118.10',
    'Bodø Lufthavn' : '118.05',
    'Bodø TMA' : '119.70',
    'Bondalen N' : '129.325',
    'Bondalen S' : '129.325',
    'Bringeland TIZ' : '118.45',
    'Brønnøy TIZ' : '119.60',
    'Båtsfjord TIZ' : '123.40',
    'Bø' : '124.35',
    'Bømoen Flyplass' : '123.50',
    'Dovre Wave 125.70' : '125.70',
    'Drangedal' : '134.05',
    'EN D102 Østre Æra' : '121.35',
    'Eggemoen A' : '120.45',
    'Eggemoen B' : '120.45',
    'Eina' : '120.45',
    'Etne' : '124.775',
    'Evenes CTR' : '119.90',
    'Evenes TMA' : '118.00',
    'Fagerhaug' : '123.50',
    'Farris TMA' : '124.35',
    'Finnmark TIA' : '126.70',
    'Flesland CTR' : '119.10',
    'Flesland TMA 2' : '125.00',
    'Flesland TMA 3' : '125.00',
    'Flesland TMA 4' : '125.00',
    'Flesland TMA 5' : '125.00',
    'Flesland TMA 6' : '121.00',
    'Flesland TMA 7' : '121.00',
    'Flesland TMA 8' : '121.00',
    'Florø TIZ' : '119.20',
    'GULLFJELLET GLIDERAREA 121.0 / 123.5' : '121.00',
    'Gardermoen CTR' : '120.10',
    'Gullknapp TIZ' : '129.900',
    'Gvarv' : '124.35',
    'Hamar Flyplass' : '130.275',
    'Hammerfest TIZ' : '121.00',
    'Hammerfest TMA' : '126.70',
    'Harstad/Narvik Lufthavn  Evenes' : '118.10',
    'Hasvik TIZ' : '119.90',
    'Helgeland TMA' : '127.90',
    'Helle TIZ' : '120.20',
    'Hokksund A' : '120.45',
    'Hokksund B' : '120.45',
    'Hoppfelt Bømoen 123.50' : '123.50',
    'Hovden TIZ' : '118.90',
    'Hvittingfoss A' : '124.35',
    'Hvittingfoss B' : '124.35',
    'Hvittingfoss C' : '124.35',
    'Hvittingfoss D' : '124.35',
    'Hvittingfoss E' : '124.35',
    'Jarlsberg' : '122.30',
    'Jotunheimen Wave 124.70' : '124.70',
    'KVAM TRANSITAREA 121.0 / 123.5' : '121.00',
    'Karmøy CTR' : '120.50',
    'Kirkenes CTR' : '120.35',
    'Kirkenes Centre TMA' : '120.35',
    'Kirkenes TMA' : '120.35',
    'Kirkenes West TMA' : '120.35',
    'Kjevik CTR' : '119.95',
    'Kjevik TMA' : '119.95',
    'Kristiansand Lufthavn Kjevik' : '119.95',
    'Kvernberget CTR' : '121.20',
    'Leknes TIZ' : '120.50',
    'Lesja Wave 125.70' : '125.70',
    'Lesja/Bjorli Flyplass' : '123.50',
    'Lofoten TMA' : '125.45',
    'Lunde' : '134.05',
    'Mehamn TIZ' : '121.20',
    'Molde Lufthavn Årø' : '119.95',
    'Molde TIZ' : '119.95',
    'Mosjøen TIZ' : '123.40',
    'Møre TMA' : '119.35',
    'Namsos TIA' : '118.550',
    'Namsos TIZ' : '119.90',
    'Notodden TIZ' : '118.80',
    'Oppdal Wave 125.70' : '125.70',
    'Oslo TMA' : '118.475',
    'Oslo TMA 2' : '120.45',
    'Oslo TMA 3' : '118.475',
    'Oslo TMA 4' : '120.45',
    'Oslo TMA 5' : '118.475',
    'Oslo TMA 6' : '118.475',
    'Oslo TMA 7' : '118.475',
    'Rena EAST 135.30' : '135.30',
    'Rena MID 135.30' : '135.30',
    'Rena WEST 135.30' : '135.30',
    'Ringebu Wave 124.775' : '124.775',
    'Rognan Flyplass' : '123.50',
    'Rondane Wave 124.775' : '124.775',
    'Rygge CTR' : '119.50',
    'Rygge Flystasjon' : '119.50',
    'Røros TIA' : '120.40',
    'Røros TIZ' : '120.40',
    'Rørvik TIZ' : '119.80',
    'Røssvoll TIZ' : '119.95',
    'Røst TIZ' : '119.25',
    'Skagen TIZ' : '120.45',
    'Skien Lufthavn  Geiteryggen' : '119.20',
    'Sogn TIA' : '124.70',
    'Sogndal TIZ' : '119.30',
    'Sola CTR' : '118.35',
    'Sola TMA' : '119.60',
    'Starmoen A' : '118.475',
    'Starmoen B' : '118.475',
    'Starmoen C' : '118.475',
    'Starmoen D' : '118.475',
    'Starmoen F' : '120.45',
    'Starmoen G' : '118.475',
    'Starmoen H' : '118.475',
    'Stokka TIZ' : '120.30',
    'Sunnvollen' : '120.45',
    'Sørkjosen TIA' : '126.70',
    'Sørkjosen TIZ' : '119.60',
    'Sørstokken TIZ' : '120.20',
    'Tokke' : '124.35',
    'Torp CTR' : '118.65',
    'Tromsø CTR' : '118.30',
    'Tromsø Lufthavn' : '118.30',
    'Tromsø TMA' : '123.75',
    'Trondheim Lufthavn  Værnes' : '122.05',
    'Tynset' : '125.70',
    'Vaagaa Wave 124.775' : '124.775',
    'Vadsø TIZ' : '118.40',
    'Valan TIZ' : '119.80',
    'Vardø TIZ' : '122.15',
    'Vigra CTR' : '119.85',
    'Værnes CTR' : '118.60',
    'Værnes TMA' : '118.60',
    'Ålesund Lufthavn Vigra' : '119.85',
    'Ørland CTR' : '118.70',
    'Ørland TMA' : '118.25',
    'Østre Æra Flyplass' : '135.30',
}

AIRSPACE_NAME = re.compile(r'((.*)(TMA|CTA|TIA|TIZ|CTR)\s+(\d*))')
def lookup_frequency(name):
    # Look up frequency in our table, from most specific to least
    if name in airspace_frequencies:   # We have frequency for exact airspace name
        return airspace_frequencies[name],name
    else:
        # For e.g. Farris TMA 3 West try in order: "Farris TMA 3", then "Farris TMA"
        if match := AIRSPACE_NAME.match(name):
            if match.group(1) in airspace_frequencies:
                return airspace_frequencies[name],match.group(1)
            else:
                name = match.group(2) + match.group(3)
                return airspace_frequencies.get(name, None),name

    return None,None
